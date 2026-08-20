/**
 *  @file   LCContent/src/LCPlugins/LCSoftwareCompensation.cc
 *
 *  @brief  Implementation of the lc software compensation plugin algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCPlugins/LCSoftwareCompensation.h"

using namespace pandora;

namespace lc_content {

LCSoftwareCompensation::LCSoftwareCompensation(const LCSoftwareCompensationParameters& lcSoftwareCompensationParameters)
    : m_softCompParameters(lcSoftwareCompensationParameters.m_softCompParameters),
      m_softCompEnergyDensityBins(lcSoftwareCompensationParameters.m_softCompEnergyDensityBins),
      m_energyDensityFinalBin(lcSoftwareCompensationParameters.m_energyDensityFinalBin),
      m_maxClusterEnergyToApplySoftComp(lcSoftwareCompensationParameters.m_maxClusterEnergyToApplySoftComp),
      m_minCleanHitEnergy(lcSoftwareCompensationParameters.m_minCleanHitEnergy),
      m_minCleanHitEnergyFraction(lcSoftwareCompensationParameters.m_minCleanHitEnergyFraction),
      m_minCleanCorrectedHitEnergy(lcSoftwareCompensationParameters.m_minCleanCorrectedHitEnergy) {
  if (9 != m_softCompParameters.size()) {
    std::cout << "LCSoftwareCompensation:LCSoftwareCompensation - Incorrect number of parameters required for software "
                 "compensation technique."
              << std::endl;
    throw STATUS_CODE_INVALID_PARAMETER;
  }
  std::sort(m_softCompEnergyDensityBins.begin(), m_softCompEnergyDensityBins.end());

  if (m_softCompEnergyDensityBins.front() < 0.f) {
    std::cout << "LCSoftwareCompensation:LCSoftwareCompensation - Input density bins contains an unphysical value"
              << std::endl;
    throw STATUS_CODE_FAILURE;
  }

  if (m_energyDensityFinalBin < m_softCompEnergyDensityBins.back()) {
    std::cout << "LCSoftwareCompensation::LCSoftwareCompensation - Input energy density final bin inconsistent with "
                 "input density bins"
              << std::endl;
    throw STATUS_CODE_FAILURE;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCSoftwareCompensation::MakeEnergyCorrections(const pandora::Cluster* const pCluster,
                                                         float& correctedHadronicEnergy) const {
  if (nullptr == pCluster || m_softCompParameters.size() != 9)
    throw STATUS_CODE_INVALID_PARAMETER;

  if (0 == pCluster->GetNCaloHits()) {
    correctedHadronicEnergy = 0.f;
    return STATUS_CODE_SUCCESS;
  }

  if (correctedHadronicEnergy < m_maxClusterEnergyToApplySoftComp) {
    const float clusterHadEnergy = pCluster->GetHadronicEnergy();
    pandora::CaloHitList clusterCaloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
    clusterCaloHitList.insert(clusterCaloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(),
                              pCluster->GetIsolatedCaloHitList().end());
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             this->SoftComp(clusterHadEnergy, clusterCaloHitList, correctedHadronicEnergy));
  }

  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CleanCluster(pCluster, correctedHadronicEnergy));

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCSoftwareCompensation::CleanCluster(const pandora::Cluster* const pCluster,
                                                float& correctedHadronicEnergy) const {
  const unsigned int firstPseudoLayer(this->GetPandora().GetPlugins()->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());
  const float clusterHadronicEnergy(pCluster->GetHadronicEnergy());

  if (std::fabs(clusterHadronicEnergy) < std::numeric_limits<float>::epsilon())
    return STATUS_CODE_FAILURE;

  bool isFineGranularity(true);
  const OrderedCaloHitList& orderedCaloHitList(pCluster->GetOrderedCaloHitList());

  for (OrderedCaloHitList::const_iterator layerIter = orderedCaloHitList.begin(),
                                          layerIterEnd = orderedCaloHitList.end();
       (layerIter != layerIterEnd) && isFineGranularity; ++layerIter) {
    const unsigned int pseudoLayer(layerIter->first);

    for (const CaloHit* const pCaloHit : *layerIter->second) {
      if (ECAL != pCaloHit->GetHitType())
        continue;

      if (this->GetPandora().GetGeometry()->GetHitTypeGranularity(pCaloHit->GetHitType()) > FINE) {
        isFineGranularity = false;
        break;
      }

      const float hitHadronicEnergy(pCaloHit->GetHadronicEnergy());

      if ((hitHadronicEnergy > m_minCleanHitEnergy) &&
          (hitHadronicEnergy / clusterHadronicEnergy > m_minCleanHitEnergyFraction)) {
        float energyInPreviousLayer(0.f);

        if (pseudoLayer > firstPseudoLayer)
          energyInPreviousLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer - 1);

        float energyInNextLayer(0.f);

        if (pseudoLayer < std::numeric_limits<unsigned int>::max())
          energyInNextLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer + 1);

        const float energyInCurrentLayer = this->GetHadronicEnergyInLayer(orderedCaloHitList, pseudoLayer);
        float energyInAdjacentLayers(energyInPreviousLayer + energyInNextLayer);

        if (pseudoLayer > firstPseudoLayer)
          energyInAdjacentLayers /= 2.f;

        float newHitHadronicEnergy(energyInAdjacentLayers - energyInCurrentLayer + hitHadronicEnergy);
        newHitHadronicEnergy = std::max(newHitHadronicEnergy, m_minCleanCorrectedHitEnergy);

        if (newHitHadronicEnergy < hitHadronicEnergy)
          correctedHadronicEnergy += newHitHadronicEnergy - hitHadronicEnergy;
      }
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LCSoftwareCompensation::GetHadronicEnergyInLayer(const OrderedCaloHitList& orderedCaloHitList,
                                                       const unsigned int pseudoLayer) const {
  OrderedCaloHitList::const_iterator iter = orderedCaloHitList.find(pseudoLayer);

  float hadronicEnergy(0.f);

  if (iter != orderedCaloHitList.end()) {
    for (const CaloHit* const pCaloHit : *iter->second) {
      hadronicEnergy += pCaloHit->GetHadronicEnergy();
    }
  }

  return hadronicEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCSoftwareCompensation::SoftComp(float clusterEnergyEstimation, const pandora::CaloHitList& caloHitList,
                                            float& energyCorrection) const {
  float energySoftComp(0.f);

  const float p1 = m_softCompParameters.at(0) + m_softCompParameters.at(1) * clusterEnergyEstimation +
                   m_softCompParameters.at(2) * clusterEnergyEstimation * clusterEnergyEstimation;
  const float p2 = m_softCompParameters.at(3) + m_softCompParameters.at(4) * clusterEnergyEstimation +
                   m_softCompParameters.at(5) * clusterEnergyEstimation * clusterEnergyEstimation;
  const float p3 = m_softCompParameters.at(6) /
                   (m_softCompParameters.at(7) + exp(m_softCompParameters.at(8) * clusterEnergyEstimation));

  for (const pandora::CaloHit* const pCaloHit : caloHitList) {
    if (HCAL == pCaloHit->GetHitType()) {
      const float hitEnergy = pCaloHit->GetHadronicEnergy();
      float rho(0.f);

      try {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindDensity(pCaloHit, rho));
      } catch (const StatusCodeException&) {
        std::cout << "LCSoftwareCompensation: Unable to find energy density of calorimeter hit" << std::endl;
      }
      const float weight(p1 * exp(p2 * rho) + p3);
      const float correctedHitEnergy(hitEnergy * weight);
      energySoftComp += correctedHitEnergy;
    } else {
      energySoftComp += pCaloHit->GetHadronicEnergy();
    }
  }

  energyCorrection = energySoftComp;

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCSoftwareCompensation::FindDensity(const pandora::CaloHit* const pCaloHit, float& energyDensity) const {
  const float mm3Todm3 = 1e-6f; // ATTN: Cell energy density defined in GeV per dm3 but Pandora cell size defined in mm,
                                // so needs conversion
  const float cellVolume =
      pCaloHit->GetCellSize0() * pCaloHit->GetCellSize1() * pCaloHit->GetCellThickness() * mm3Todm3;

  if (cellVolume < std::numeric_limits<float>::epsilon()) {
    throw STATUS_CODE_FAILURE;
  }

  const float hitEnergyHadronic(pCaloHit->GetHadronicEnergy());
  const float hitEnergyDensity(hitEnergyHadronic / cellVolume);

  if (hitEnergyDensity >= m_softCompEnergyDensityBins.back()) {
    energyDensity = m_energyDensityFinalBin;
    return STATUS_CODE_SUCCESS;
  } else {
    for (unsigned int iBin = 0; iBin < m_softCompEnergyDensityBins.size() - 1; iBin++) {
      const float lowBinContent = m_softCompEnergyDensityBins.at(iBin);
      const float highBinContent = m_softCompEnergyDensityBins.at(iBin + 1);

      if (hitEnergyDensity >= lowBinContent && hitEnergyDensity < highBinContent) {
        energyDensity = (lowBinContent + highBinContent) * 0.5f;
        return STATUS_CODE_SUCCESS;
      }
    }
  }

  std::cout << "LCSoftwareCompensation::FindDensity - EnergyDensity binning inconsistency" << std::endl;

  return STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCSoftwareCompensation::ReadSettings(const TiXmlHandle /*xmlHandle*/) { return STATUS_CODE_SUCCESS; }

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LCSoftwareCompensationParameters::LCSoftwareCompensationParameters()
    : m_softCompParameters(), m_softCompEnergyDensityBins(), m_energyDensityFinalBin(30.f),
      m_maxClusterEnergyToApplySoftComp(100.f), m_minCleanHitEnergy(0.5f), m_minCleanHitEnergyFraction(0.01f),
      m_minCleanCorrectedHitEnergy(0.1f) {
  const unsigned int nParameters(9);
  const float parameters[nParameters] = {2.49632f,      -0.0697302f, 0.000946986f, -0.112311f, 0.0028182f,
                                         -9.62602e-05f, 0.168614f,   0.224318f,    -0.0872853f};
  m_softCompParameters.insert(m_softCompParameters.begin(), parameters, parameters + nParameters);

  const unsigned int nBins(10);
  const float bins[nBins] = {0.f, 2.f, 5.f, 7.5f, 9.5f, 13.f, 16.f, 20.f, 23.5f, 28.f};
  m_softCompEnergyDensityBins.insert(m_softCompEnergyDensityBins.begin(), bins, bins + nBins);
}

} // namespace lc_content
