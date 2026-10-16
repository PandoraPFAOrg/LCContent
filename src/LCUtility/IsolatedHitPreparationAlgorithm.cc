#include "LCUtility/IsolatedHitPreparationAlgorithm.h"

#include "Api/PandoraContentApi.h"
#include "Helpers/XmlHelper.h"
#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"

#include <unordered_set>

using namespace pandora;

namespace lc_content {

Algorithm* IsolatedHitPreparationAlgorithm::Factory::CreateAlgorithm() const {
  return new IsolatedHitPreparationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitPreparationAlgorithm::Run() {
  // retrieve clusters
  const pandora::ClusterList* pClusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
  // retrieve calo hits
  const pandora::CaloHitList* pCaloHitList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

  // collect hits that are not clustered
  std::unordered_set<const pandora::CaloHit*> clusteredHits;
  for (const auto* pCluster : *pClusterList) {
    pandora::CaloHitList clusterHits;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHits);
    for (const auto* pHit : clusterHits)
      clusteredHits.insert(pHit);
  } // loop over clusters

  // set hits that are not clustered as isolated
  for (const auto* pHit : *pCaloHitList) {
    if (clusteredHits.count(pHit) == 0) {
      if (m_excludeEcal && pHit->GetElectromagneticEnergy() > 0.f)
        continue; // configured to exclude ECAL hits from isolated hit flagging
      if (m_excludeHcal && pHit->GetHadronicEnergy() > 0.f)
        continue; // configured to exclude HCAL hits from isolated hit flagging

      // energy zero-suppression: drop isolated-hit candidates below the configured cut.
      // The cut is a direct energy (GeV) threshold -- set it to the energy that corresponds
      // to your desired photon count.  A separate cherenkov cut is applied only when
      // CherEnergyThreshold > 0; otherwise the (scintillation) cut applies to all.
      const bool isCher = pHit->GetHitType() == pandora::DRC_CHEREN;
      const float threshold = (isCher && m_cherEnergyThreshold > 0.f) ? m_cherEnergyThreshold : m_scintEnergyThreshold;
      if (pHit->GetInputEnergy() <= threshold)
        continue; // sub-threshold -> never flagged isolated -> not merged / not regressed

      PandoraContentApi::CaloHit::Metadata metadata;
      metadata.m_isIsolated = true;
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::CaloHit::AlterMetadata(*this, pHit, metadata));
    }
  } // loop over calo hits

  return pandora::STATUS_CODE_SUCCESS;
} // Run

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "ExcludeEcalHits", m_excludeEcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "ExcludeHcalHits", m_excludeHcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "ScintEnergyThreshold", m_scintEnergyThreshold));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "CherEnergyThreshold", m_cherEnergyThreshold));

  return pandora::STATUS_CODE_SUCCESS;
} // ReadSettings

} // namespace lc_content
