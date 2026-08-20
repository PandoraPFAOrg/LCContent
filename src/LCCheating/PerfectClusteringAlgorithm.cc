/**
 *  @file   LCContent/src/LCCheating/PerfectClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the cheating clustering algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCCheating/PerfectClusteringAlgorithm.h"

using namespace pandora;

namespace lc_content {

PerfectClusteringAlgorithm::PerfectClusteringAlgorithm()
    : m_shouldUseOnlyECalHits(false), m_shouldUseIsolatedHits(false), m_simpleMCParticleCollection(true),
      m_minWeightFraction(0.01f) {}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PerfectClusteringAlgorithm::SelectMCParticlesForClustering(const MCParticle* const pMCParticle) const {
  if (m_particleIdList.empty())
    return true;

  for (IntVector::const_iterator iter = m_particleIdList.begin(), iterEnd = m_particleIdList.end(); iter != iterEnd;
       ++iter) {
    if (pMCParticle->GetParticleId() == *iter)
      return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PerfectClusteringAlgorithm::Run() {
  const CaloHitList* pCaloHitList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

  CaloHitList localCaloHitList(pCaloHitList->begin(), pCaloHitList->end());
  MCParticleToHitListMap mcParticleToHitListMap;

  for (CaloHitList::const_iterator hitIter = localCaloHitList.begin(), hitIterEnd = localCaloHitList.end();
       hitIter != hitIterEnd; ++hitIter) {
    try {
      const CaloHit* const pCaloHit = *hitIter;

      if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
        continue;

      if (m_shouldUseOnlyECalHits && (ECAL != pCaloHit->GetHitType()))
        continue;

      if (!m_shouldUseIsolatedHits && pCaloHit->IsIsolated())
        continue;

      if (m_simpleMCParticleCollection) {
        this->SimpleMCParticleCollection(pCaloHit, mcParticleToHitListMap);
      } else {
        this->FullMCParticleCollection(pCaloHit, mcParticleToHitListMap);
      }
    } catch (StatusCodeException&) {
    }
  }

  this->CreateClusters(mcParticleToHitListMap);

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PerfectClusteringAlgorithm::SimpleMCParticleCollection(const CaloHit* const pCaloHit,
                                                            MCParticleToHitListMap& mcParticleToHitListMap) const {
  const MCParticle* const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

  if (!this->SelectMCParticlesForClustering(pMCParticle))
    return;

  this->AddToHitListMap(pCaloHit, pMCParticle, mcParticleToHitListMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PerfectClusteringAlgorithm::FullMCParticleCollection(const CaloHit* const pCaloHit,
                                                          MCParticleToHitListMap& mcParticleToHitListMap) const {
  const MCParticleWeightMap mcParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

  if (mcParticleWeightMap.size() < 2)
    return this->SimpleMCParticleCollection(pCaloHit, mcParticleToHitListMap);

  float mcParticleWeightSum(0.f);

  MCParticleList mcParticleList;
  for (const auto& mapEntry : mcParticleWeightMap)
    mcParticleList.push_back(mapEntry.first);
  mcParticleList.sort(PointerLessThan<MCParticle>());

  for (const MCParticle* const pMCParticle : mcParticleList)
    mcParticleWeightSum += mcParticleWeightMap.at(pMCParticle);

  if (mcParticleWeightSum < std::numeric_limits<float>::epsilon())
    throw StatusCodeException(STATUS_CODE_FAILURE);

  const CaloHit* pLocalCaloHit = pCaloHit;

  for (const MCParticle* const pMCParticle : mcParticleList) {
    const float weight(mcParticleWeightMap.at(pMCParticle));

    if (!this->SelectMCParticlesForClustering(pMCParticle))
      continue;

    const CaloHit* pCaloHitToAdd = pLocalCaloHit;

    if (pCaloHitToAdd->GetWeight() < std::numeric_limits<float>::epsilon())
      throw StatusCodeException(STATUS_CODE_FAILURE);

    const float weightFraction(weight / (mcParticleWeightSum * pCaloHitToAdd->GetWeight()));
    const bool isBelowThreshold((weightFraction - m_minWeightFraction) < std::numeric_limits<float>::epsilon());

    if (isBelowThreshold)
      continue;

    const bool shouldFragment(weightFraction + m_minWeightFraction - 1.f < std::numeric_limits<float>::epsilon());

    if (shouldFragment)
      PANDORA_THROW_RESULT_IF(
          STATUS_CODE_SUCCESS, !=,
          PandoraContentApi::Fragment(*this, pLocalCaloHit, weightFraction, pCaloHitToAdd, pLocalCaloHit));

    this->AddToHitListMap(pCaloHitToAdd, pMCParticle, mcParticleToHitListMap);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PerfectClusteringAlgorithm::AddToHitListMap(const CaloHit* const pCaloHitToAdd,
                                                 const MCParticle* const pMCParticle,
                                                 MCParticleToHitListMap& mcParticleToHitListMap) const {
  MCParticleToHitListMap::iterator iter(mcParticleToHitListMap.find(pMCParticle));

  if (mcParticleToHitListMap.end() == iter) {
    CaloHitList* const pCaloHitList = new CaloHitList(1, pCaloHitToAdd);
    (void)mcParticleToHitListMap.insert(MCParticleToHitListMap::value_type(pMCParticle, pCaloHitList));
  } else {
    iter->second->push_back(pCaloHitToAdd);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PerfectClusteringAlgorithm::CreateClusters(const MCParticleToHitListMap& mcParticleToHitListMap) const {
  MCParticleList mcParticleList;
  for (const auto& mapEntry : mcParticleToHitListMap)
    mcParticleList.push_back(mapEntry.first);
  mcParticleList.sort(PointerLessThan<MCParticle>());

  for (const MCParticle* const pMCParticle : mcParticleList) {
    CaloHitList* const pCaloHitList(mcParticleToHitListMap.at(pMCParticle));

    if (!pCaloHitList->empty()) {
      const Cluster* pCluster = NULL;
      PandoraContentApi::Cluster::Parameters parameters;
      parameters.m_caloHitList = *pCaloHitList;
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

      PandoraContentApi::Cluster::Metadata metadata;

      switch (pMCParticle->GetParticleId()) {
      case PHOTON:
      case E_PLUS:
      case E_MINUS:
      case MU_PLUS:
      case MU_MINUS:
        metadata.m_particleId = pMCParticle->GetParticleId();
        break;
      default:
        break;
      }

      if (metadata.m_particleId.IsInitialized())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
    }
    delete pCaloHitList;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PerfectClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadVectorOfValues(xmlHandle, "ParticleIdList", m_particleIdList));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ShouldUseOnlyECalHits", m_shouldUseOnlyECalHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ShouldUseIsolatedHits", m_shouldUseIsolatedHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "SimpleMCParticleCollection", m_simpleMCParticleCollection));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinWeightFraction", m_minWeightFraction));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
