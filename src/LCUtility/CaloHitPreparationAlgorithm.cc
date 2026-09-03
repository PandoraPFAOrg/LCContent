/**
 *  @file   LCContent/src/LCUtility/CaloHitPreparationAlgorithm.cc
 *
 *  @brief  Implementation of the calo hit preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCUtility/CaloHitPreparationAlgorithm.h"
#include "LCUtility/KDTreeLinkerAlgoT.h"

#include <fstream>
#include <sstream>

using namespace pandora;

namespace lc_content {

CaloHitPreparationAlgorithm::CaloHitPreparationAlgorithm()
    : m_caloHitMaxSeparation2(100.f * 100.f), m_isolationCaloHitMaxSeparation2(1000.f * 1000.f), m_isolationNLayers(2),
      m_isolationCutDistanceFine2(25.f * 25.f), m_isolationCutDistanceCoarse2(200.f * 200.f),
      m_isolationSearchSafetyFactor(2.f), m_isolationMaxNearbyHits(2), m_mipLikeMipCut(5.f), m_mipNCellsForNearbyHit(2),
      m_mipMaxNearbyHits(1), m_bibNLayers(2), m_bibCaloHitMaxSeparationSquared(1000.f * 1000.f),
      m_bibCutDistanceSquared(200.f * 200.f), m_bibSearchSafetyFactor(2.f), m_bibTimingCutEnabled(false),
      m_bibTimingWindowLow(0.f), m_bibTimingWindowHigh(0.f), m_hitNodes4D(new std::vector<HitKDNode4D>),
      m_hitsKdTree4D(new HitKDTree4D) {}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitPreparationAlgorithm::~CaloHitPreparationAlgorithm() {
  delete m_hitNodes4D;
  delete m_hitsKdTree4D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitPreparationAlgorithm::Run() {
  try {
    const CaloHitList* pCaloHitList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    this->InitializeKDTree(pCaloHitList);

    OrderedCaloHitList orderedCaloHitList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end();
         iter != iterEnd; ++iter) {
      for (CaloHitList::iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end();
           hitIter != hitIterEnd; ++hitIter) {
        this->CalculateCaloHitProperties(*hitIter, orderedCaloHitList);
      }
    }
  } catch (StatusCodeException& statusCodeException) {
    std::cout << "CaloHitPreparationAlgorithm: Failed to calculate calo hit properties, "
              << statusCodeException.ToString() << std::endl;
    return statusCodeException.GetStatusCode();
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitPreparationAlgorithm::InitializeKDTree(const CaloHitList* const pCaloHitList) {
  m_hitsKdTree4D->clear();
  m_hitNodes4D->clear();
  KDTreeTesseract hitsBoundingRegion4D = fill_and_bound_4d_kd_tree(this, *pCaloHitList, *m_hitNodes4D, true);
  m_hitsKdTree4D->build(*m_hitNodes4D, hitsBoundingRegion4D);
  m_hitNodes4D->clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitPreparationAlgorithm::CalculateCaloHitProperties(const CaloHit* const pCaloHit,
                                                             const OrderedCaloHitList& orderedCaloHitList) {
  // Calculate number of adjacent pseudolayers to examine
  const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());
  const unsigned int isolationMaxLayer(pseudoLayer + m_isolationNLayers);
  const unsigned int isolationMinLayer((pseudoLayer < m_isolationNLayers) ? 0 : pseudoLayer - m_isolationNLayers);

  // Initialize variables
  bool isIsolated = true;
  unsigned int isolationNearbyHits = 0;

  // Loop over adjacent pseudolayers
  for (unsigned int iPseudoLayer = isolationMinLayer; iPseudoLayer <= isolationMaxLayer; ++iPseudoLayer) {
    OrderedCaloHitList::const_iterator adjacentPseudoLayerIter = orderedCaloHitList.find(iPseudoLayer);

    if (orderedCaloHitList.end() == adjacentPseudoLayerIter)
      continue;

    // IsIsolated flag
    if (isIsolated && (isolationMinLayer <= iPseudoLayer) && (isolationMaxLayer >= iPseudoLayer)) {
      isolationNearbyHits += this->IsolationCountNearbyHits(iPseudoLayer, pCaloHit);
      isIsolated = isolationNearbyHits < m_isolationMaxNearbyHits;
    }

    // Possible mip flag
    if (pseudoLayer == iPseudoLayer) {
      if (MUON == pCaloHit->GetHitType()) {
        PandoraContentApi::CaloHit::Metadata metadata;
        metadata.m_isPossibleMip = true;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
        continue;
      }

      const CartesianVector& positionVector(pCaloHit->GetPositionVector());

      const float x(positionVector.GetX());
      const float y(positionVector.GetY());

      const float angularCorrection((BARREL == pCaloHit->GetHitRegion())
                                        ? positionVector.GetMagnitude() / std::sqrt(x * x + y * y)
                                        : positionVector.GetMagnitude() / std::fabs(positionVector.GetZ()));

      if ((pCaloHit->GetMipEquivalentEnergy() <= (m_mipLikeMipCut * angularCorrection) || pCaloHit->IsDigital()) &&
          (m_mipMaxNearbyHits >= this->MipCountNearbyHits(iPseudoLayer, pCaloHit))) {
        PandoraContentApi::CaloHit::Metadata metadata;
        metadata.m_isPossibleMip = true;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
      }
    }
  }

  if (isIsolated) {
    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = true;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
  }

  // IsPossibleBIB flag: nearby energy density below the per-pseudolayer cut loaded for this
  // hit's subdetector and region, and/or (when enabled) a hit time outside the accepted window.
  // A zero cut can never flag, so the neighbourhood sum is skipped in that case.
  bool isPossibleBIB = false;

  const float bibEnergyDensityCut(this->GetBIBEnergyDensityCut(pCaloHit));

  if (bibEnergyDensityCut > 0.f) {
    const unsigned int bibMaxLayer(pseudoLayer + m_bibNLayers);
    const unsigned int bibMinLayer((pseudoLayer < m_bibNLayers) ? 0 : pseudoLayer - m_bibNLayers);

    float bibNearbyEnergy = 0.f;

    for (unsigned int iPseudoLayer = bibMinLayer; iPseudoLayer <= bibMaxLayer; ++iPseudoLayer) {
      if (orderedCaloHitList.end() == orderedCaloHitList.find(iPseudoLayer))
        continue;

      bibNearbyEnergy += this->BIBSumNearbyEnergy(iPseudoLayer, pCaloHit);
    }

    // Density in GeV/mm^2/layer: summed energy over the sampled neighbourhood, normalized to its
    // transverse area and the number of pseudolayers examined
    const unsigned int bibNSampledLayers(bibMaxLayer - bibMinLayer + 1);
    const float bibEnergyDensity(
        bibNearbyEnergy / (static_cast<float>(M_PI) * m_bibCutDistanceSquared * static_cast<float>(bibNSampledLayers)));

    isPossibleBIB = (bibEnergyDensity < bibEnergyDensityCut);
  }

  if (m_bibTimingCutEnabled &&
      ((pCaloHit->GetTime() < m_bibTimingWindowLow) || (pCaloHit->GetTime() > m_bibTimingWindowHigh)))
    isPossibleBIB = true;

  if (isPossibleBIB) {
    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isPossibleBIB = true;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                            PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int CaloHitPreparationAlgorithm::IsolationCountNearbyHits(unsigned int searchLayer,
                                                                   const CaloHit* const pCaloHit) {
  const CartesianVector& positionVector(pCaloHit->GetPositionVector());
  const float positionMagnitudeSquared(positionVector.GetMagnitudeSquared());
  const float isolationCutDistanceSquared(
      (PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pCaloHit->GetHitType()) <= FINE)
          ? m_isolationCutDistanceFine2
          : m_isolationCutDistanceCoarse2);

  unsigned int nearbyHitsFound = 0;

  // construct the kd tree search
  CaloHitList nearby_hits;
  const float searchDistance(m_isolationSearchSafetyFactor * std::sqrt(isolationCutDistanceSquared));
  KDTreeTesseract searchRegionHits =
      build_4d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance, searchLayer);

  std::vector<HitKDNode4D> found;
  m_hitsKdTree4D->search(searchRegionHits, found);

  for (const auto& hit : found) {
    nearby_hits.push_back(hit.data);
  }

  for (CaloHitList::const_iterator iter = nearby_hits.begin(), iterEnd = nearby_hits.end(); iter != iterEnd; ++iter) {
    if (pCaloHit == *iter)
      continue;

    const CartesianVector positionDifference(positionVector - (*iter)->GetPositionVector());
    const CartesianVector crossProduct(positionVector.GetCrossProduct(positionDifference));

    if (positionDifference.GetMagnitudeSquared() > m_isolationCaloHitMaxSeparation2)
      continue;

    if ((crossProduct.GetMagnitudeSquared() / positionMagnitudeSquared) < isolationCutDistanceSquared)
      ++nearbyHitsFound;
  }

  return nearbyHitsFound;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int CaloHitPreparationAlgorithm::MipCountNearbyHits(unsigned int searchLayer, const CaloHit* const pCaloHit) {
  const float mipNCellsForNearbyHit(m_mipNCellsForNearbyHit + 0.5f);

  unsigned int nearbyHitsFound = 0;
  const CartesianVector& positionVector(pCaloHit->GetPositionVector());
  const bool isHitInBarrelRegion(pCaloHit->GetHitRegion() == BARREL);

  // construct the kd tree search
  CaloHitList nearby_hits;
  const float searchDistance(std::sqrt(m_caloHitMaxSeparation2));
  KDTreeTesseract searchRegionHits =
      build_4d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance, searchLayer);

  std::vector<HitKDNode4D> found;
  m_hitsKdTree4D->search(searchRegionHits, found);

  for (const auto& hit : found) {
    nearby_hits.push_back(hit.data);
  }

  for (CaloHitList::const_iterator iter = nearby_hits.begin(), iterEnd = nearby_hits.end(); iter != iterEnd; ++iter) {
    if (pCaloHit == *iter)
      continue;

    const CartesianVector positionDifference(positionVector - (*iter)->GetPositionVector());

    if (positionDifference.GetMagnitudeSquared() > m_caloHitMaxSeparation2)
      continue;

    const float cellLengthScale(pCaloHit->GetCellLengthScale());

    if (isHitInBarrelRegion) {
      const float dX(std::fabs(positionDifference.GetX()));
      const float dY(std::fabs(positionDifference.GetY()));
      const float dZ(std::fabs(positionDifference.GetZ()));
      const float dPhi(std::sqrt(dX * dX + dY * dY));

      if ((dZ < (mipNCellsForNearbyHit * cellLengthScale)) && (dPhi < (mipNCellsForNearbyHit * cellLengthScale)))
        ++nearbyHitsFound;
    } else {
      const float dX(std::fabs(positionDifference.GetX()));
      const float dY(std::fabs(positionDifference.GetY()));

      if ((dX < (mipNCellsForNearbyHit * cellLengthScale)) && (dY < (mipNCellsForNearbyHit * cellLengthScale)))
        ++nearbyHitsFound;
    }
  }

  return nearbyHitsFound;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloHitPreparationAlgorithm::BIBSumNearbyEnergy(unsigned int searchLayer, const CaloHit* const pCaloHit) {
  const CartesianVector& positionVector(pCaloHit->GetPositionVector());
  const float positionMagnitudeSquared(positionVector.GetMagnitudeSquared());

  float nearbyEnergy = 0.f;

  // construct the kd tree search
  CaloHitList nearby_hits;
  const float searchDistance(m_bibSearchSafetyFactor * std::sqrt(m_bibCutDistanceSquared));
  KDTreeTesseract searchRegionHits =
      build_4d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance, searchLayer);

  std::vector<HitKDNode4D> found;
  m_hitsKdTree4D->search(searchRegionHits, found);

  for (const auto& hit : found) {
    nearby_hits.push_back(hit.data);
  }

  for (CaloHitList::const_iterator iter = nearby_hits.begin(), iterEnd = nearby_hits.end(); iter != iterEnd; ++iter) {
    // The hit's own energy is deliberately excluded: including it lets the most
    // energetic BIB hits lift themselves above the cut, notably weakening the veto
    if (pCaloHit == *iter)
      continue;

    const CartesianVector positionDifference(positionVector - (*iter)->GetPositionVector());
    const CartesianVector crossProduct(positionVector.GetCrossProduct(positionDifference));

    if (positionDifference.GetMagnitudeSquared() > m_bibCaloHitMaxSeparationSquared)
      continue;

    if ((crossProduct.GetMagnitudeSquared() / positionMagnitudeSquared) < m_bibCutDistanceSquared)
      nearbyEnergy += (*iter)->GetElectromagneticEnergy();
  }

  return nearbyEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CaloHitPreparationAlgorithm::GetBIBEnergyDensityCut(const CaloHit* const pCaloHit) const {
  const bool isEndcap(ENDCAP == pCaloHit->GetHitRegion());
  const FloatVector* pCuts(NULL);

  switch (pCaloHit->GetHitType()) {
  case ECAL:
    pCuts = isEndcap ? &m_bibCutsECalEndcap : &m_bibCutsECalBarrel;
    break;
  case HCAL:
    pCuts = isEndcap ? &m_bibCutsHCalEndcap : &m_bibCutsHCalBarrel;
    break;
  default:
    return 0.f;
  }

  const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());
  return (pseudoLayer < pCuts->size()) ? (*pCuts)[pseudoLayer] : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitPreparationAlgorithm::LoadBIBCuts(const std::string& fileName, FloatVector& cuts) const {
  if (fileName.empty())
    return STATUS_CODE_SUCCESS;

  std::ifstream file(fileName.c_str());

  if (!file.is_open()) {
    std::cout << "CaloHitPreparationAlgorithm: could not open BIB cuts file " << fileName << std::endl;
    return STATUS_CODE_NOT_FOUND;
  }

  std::string line;

  while (std::getline(file, line)) {
    const std::size_t comment(line.find('#'));

    if (std::string::npos != comment)
      line.resize(comment);

    std::istringstream parser(line);
    unsigned int pseudoLayer(0);
    float cut(0.f);

    if (!(parser >> pseudoLayer >> cut))
      continue;

    if (cuts.size() <= pseudoLayer)
      cuts.resize(pseudoLayer + 1, 0.f);

    cuts[pseudoLayer] = cut;
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  float caloHitMaxSeparation(std::sqrt(m_caloHitMaxSeparation2));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "CaloHitMaxSeparation", caloHitMaxSeparation));
  m_caloHitMaxSeparation2 = caloHitMaxSeparation * caloHitMaxSeparation;

  float isolationCaloHitMaxSeparation(std::sqrt(m_isolationCaloHitMaxSeparation2));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "IsolationCaloHitMaxSeparation", isolationCaloHitMaxSeparation));
  m_isolationCaloHitMaxSeparation2 = isolationCaloHitMaxSeparation * isolationCaloHitMaxSeparation;

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "IsolationNLayers", m_isolationNLayers));

  float isolationCutDistanceFine(std::sqrt(m_isolationCutDistanceFine2));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "IsolationCutDistanceFine", isolationCutDistanceFine));
  m_isolationCutDistanceFine2 = isolationCutDistanceFine * isolationCutDistanceFine;

  float isolationCutDistanceCoarse(std::sqrt(m_isolationCutDistanceCoarse2));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "IsolationCutDistanceCoarse", isolationCutDistanceCoarse));
  m_isolationCutDistanceCoarse2 = isolationCutDistanceCoarse * isolationCutDistanceCoarse;

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "IsolationSearchSafetyFactor", m_isolationSearchSafetyFactor));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "IsolationMaxNearbyHits", m_isolationMaxNearbyHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MipLikeMipCut", m_mipLikeMipCut));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MipNCellsForNearbyHit", m_mipNCellsForNearbyHit));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MipMaxNearbyHits", m_mipMaxNearbyHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBNLayers", m_bibNLayers));

  float bibCaloHitMaxSeparation(std::sqrt(m_bibCaloHitMaxSeparationSquared));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCaloHitMaxSeparation", bibCaloHitMaxSeparation));
  m_bibCaloHitMaxSeparationSquared = bibCaloHitMaxSeparation * bibCaloHitMaxSeparation;

  float bibCutDistance(std::sqrt(m_bibCutDistanceSquared));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCutDistance", bibCutDistance));
  m_bibCutDistanceSquared = bibCutDistance * bibCutDistance;

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBSearchSafetyFactor", m_bibSearchSafetyFactor));

  // Per-pseudolayer energy density cuts, one optional file per subdetector and region; a
  // subdetector/region without a file is disabled (no hits there are ever flagged)
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCutsFileECalBarrel", m_bibCutsFileECalBarrel));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCutsFileECalEndcap", m_bibCutsFileECalEndcap));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCutsFileHCalBarrel", m_bibCutsFileHCalBarrel));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBCutsFileHCalEndcap", m_bibCutsFileHCalEndcap));

  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->LoadBIBCuts(m_bibCutsFileECalBarrel, m_bibCutsECalBarrel));
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->LoadBIBCuts(m_bibCutsFileECalEndcap, m_bibCutsECalEndcap));
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->LoadBIBCuts(m_bibCutsFileHCalBarrel, m_bibCutsHCalBarrel));
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->LoadBIBCuts(m_bibCutsFileHCalEndcap, m_bibCutsHCalEndcap));

  // Optional timing window: when enabled, a hit with time outside [low, high] ns is flagged
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBTimingCutEnabled", m_bibTimingCutEnabled));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBTimingWindowLow", m_bibTimingWindowLow));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BIBTimingWindowHigh", m_bibTimingWindowHigh));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
