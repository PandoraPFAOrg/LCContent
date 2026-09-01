/**
 *  @file   PandoraSDK/src/Templates/CheatingParticleIDAlgorithm.cc
 *
 *  @brief  Implementation of the cheating particle id algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCCheating/CheatingParticleIDAlgorithm.h"

using namespace pandora;

CheatingParticleIDAlgorithm::CheatingParticleIDAlgorithm() : m_useClusterOverTrackID(true) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingParticleIDAlgorithm::Run() {
  const PfoList* pPfoList(nullptr);
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pPfoList));

  for (const ParticleFlowObject* const pPfo : *pPfoList) {
    try {
      const ClusterList clusterList = pPfo->GetClusterList();
      const TrackList trackList = pPfo->GetTrackList();

      const int nClusters(clusterList.size());
      const int nTracks(trackList.size());

      if (nTracks > 1) {
        std::cout << "Unable to associate a single mc particle to a pfo with multiple tracks." << std::endl;
        continue;
      }

      PandoraContentApi::ParticleFlowObject::Metadata metadata;

      // Case 1) No Clusters or Tracks
      if (nClusters == 0 && nTracks == 0) {
        throw StatusCodeException(STATUS_CODE_FAILURE);
      }
      // Case 2) One Track Only
      else if (nClusters == 0 && nTracks == 1) {
        const Track* const pTrack(trackList.front());
        const MCParticle* const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pTrack));
        metadata.m_particleId = pMainMCParticle->GetParticleId();
      }
      // Case 3) Clusters With No Tracks
      else if (nClusters > 0 && nTracks == 0) {
        const MCParticle* const pMainMCParticle(MCParticleHelper::GetMainMCParticle(&clusterList));
        metadata.m_particleId = pMainMCParticle->GetParticleId();
      }
      // Case 4) Clusters and Tracks
      else if (nClusters > 0 && nTracks == 1) {
        const MCParticle* const pMainClustersMCParticle(MCParticleHelper::GetMainMCParticle(&clusterList));
        const int particleIDCluster(pMainClustersMCParticle->GetParticleId());

        const Track* const pTrack(trackList.front());
        const MCParticle* const pMainTrackMCParticle(MCParticleHelper::GetMainMCParticle(pTrack));
        const int particleIDTrack(pMainTrackMCParticle->GetParticleId());

        if (particleIDCluster == particleIDTrack) {
          metadata.m_particleId = particleIDCluster;
        } else if (m_useClusterOverTrackID) {
          std::cout << "Conflict between track and cluster best mc particle.  Using cluster mc particle." << std::endl;
          metadata.m_particleId = particleIDCluster;
        } else {
          std::cout << "Conflict between track and cluster best mc particle.  Using track mc particle." << std::endl;
          metadata.m_particleId = particleIDTrack;
        }
      }

      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
    } catch (StatusCodeException& statusCodeException) {
      std::cout << "CheatingParticleIDAlgorithm::Run unable to associate mc particle to pfo: "
                << statusCodeException.ToString() << std::endl;
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingParticleIDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "UseClusterOverTrackID", m_useClusterOverTrackID));

  return STATUS_CODE_SUCCESS;
}
