/**
 *  @file   LCContent/src/LCPlugins/DualReadoutCorrection.cc
 *
 *  @brief  Implementation of the dual-readout correction plugin algorithm class.
 *
 *  $Log: $
 */

#include "LCPlugins/DualReadoutCorrection.h"
#include "Pandora/AlgorithmHeaders.h"

namespace lc_content {

DualReadoutEnergy RunDualReadoutCorrection(const pandora::Cluster* const aClus, const float chiEcal,
                                           const float chiHcal) {
  // initialize output struct
  DualReadoutEnergy out;

  pandora::CaloHitList pandoraCaloHitList;
  aClus->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
  const auto& isolatedHits = aClus->GetIsolatedCaloHitList();

  if (pandoraCaloHitList.empty())
    return out;

  for (const pandora::CaloHitList& hitList : {pandoraCaloHitList, isolatedHits}) {
    for (const auto* hit : hitList) {
      float ecalEnergy = hit->GetElectromagneticEnergy();
      float hcalEnergy = hit->GetHadronicEnergy();

      if (ecalEnergy > 0.f) {
        if (hit->GetHitType() == pandora::DRC_SCINT) {
          out.ecalS += ecalEnergy;
          out.ecalS_byDepth[hit->GetLayer()] += ecalEnergy;
        } else if (hit->GetHitType() == pandora::DRC_CHEREN) {
          out.ecalC += ecalEnergy;
          out.ecalC_byDepth[hit->GetLayer()] += ecalEnergy;
        }
      }

      if (hcalEnergy > 0.f) {
        if (hit->GetHitType() == pandora::DRC_SCINT)
          out.hcalS += hcalEnergy;
        else if (hit->GetHitType() == pandora::DRC_CHEREN)
          out.hcalC += hcalEnergy;
      }
    } // hit loop
  } // loop over hit lists (regular hits and isolated hits)

  // apply dual-readout correction
  // Note: we assume that chiEcal and chiHcal inputs are checked to be less than 1
  // so that the denominators here are not close to zero
  // (it should be checked in the ReadSettings of the plugin)
  out.energyEcalCorrected = (out.ecalS - chiEcal * out.ecalC) / (1.f - chiEcal);
  out.energyHcalCorrected = (out.hcalS - chiHcal * out.hcalC) / (1.f - chiHcal);

  return out;
} // RunDualReadoutCorrection

pandora::StatusCode DualReadoutCorrection::MakeEnergyCorrections(const pandora::Cluster* const pCluster,
                                                                 float& correctedEnergy) const {
  auto dualReadoutEnergy(RunDualReadoutCorrection(pCluster, m_chiEcal, m_chiHcal));

  // use the sum of corrected ECAL and HCAL energy as the corrected energy
  correctedEnergy = dualReadoutEnergy.energyEcalCorrected + dualReadoutEnergy.energyHcalCorrected;

  return pandora::STATUS_CODE_SUCCESS;
}

pandora::StatusCode DualReadoutCorrection::ReadSettings(const pandora::TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                           pandora::XmlHelper::ReadValue(xmlHandle, "ChiEcal", m_chiEcal));

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                           pandora::XmlHelper::ReadValue(xmlHandle, "ChiHcal", m_chiHcal));

  if (m_chiEcal >= 1.f || m_chiHcal >= 1.f || m_chiEcal <= 0.f || m_chiHcal <= 0.f) {
    std::cerr << "DualReadoutCorrection::ReadSettings: ChiEcal and ChiHcal must be 0 < chi < 1";
    std::cerr << " to avoid singularity in dual-readout correction." << std::endl;

    return pandora::STATUS_CODE_INVALID_PARAMETER;
  }

  return pandora::STATUS_CODE_SUCCESS;
} // ReadSettings

} // namespace lc_content