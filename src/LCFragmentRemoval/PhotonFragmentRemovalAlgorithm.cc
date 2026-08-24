/**
 *  @file   LCContent/src/LCFragmentRemoval/PhotonFragmentRemovalAlgorithm.cc
 *
 *  @brief  Implementation of the photon fragment removal algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCFragmentRemoval/PhotonFragmentRemovalAlgorithm.h"

#include "LCHelpers/SortingHelper.h"

using namespace pandora;

namespace lc_content {

PhotonFragmentRemovalAlgorithm::PhotonFragmentRemovalAlgorithm()
    : m_nMaxPasses(200), m_innerLayerTolerance(5), m_minCosOpeningAngle(0.95f), m_useOnlyPhotonLikeDaughters(true),
      m_photonLikeMaxInnerLayer(10), m_photonLikeMinDCosR(0.5f), m_photonLikeMaxShowerStart(5.f),
      m_photonLikeMaxProfileDiscrepancy(0.75f), m_contactCutNLayers(2), m_contactCutConeFraction1(0.5f),
      m_contactCutCloseHitFraction1(0.5f), m_contactCutCloseHitFraction2(0.2f), m_contactEvidenceNLayers(2),
      m_contactEvidenceFraction(0.5f), m_coneEvidenceFraction1(0.5f), m_distanceEvidence1(100.f),
      m_distanceEvidence1d(100.f), m_distanceEvidenceCloseFraction1Multiplier(1.f),
      m_distanceEvidenceCloseFraction2Multiplier(2.f), m_contactWeight(1.f), m_coneWeight(1.f), m_distanceWeight(1.f),
      m_minEvidence(2.f) {
  m_minDaughterCaloHits = 5;
  m_minDaughterHadronicEnergy = 0.025f;
  m_contactCutMaxDistance = 20.f;

  m_contactParameters.m_coneCosineHalfAngle1 = 0.95f;
  m_contactParameters.m_closeHitDistance1 = 40.f;
  m_contactParameters.m_closeHitDistance2 = 20.f;
  m_contactParameters.m_minCosOpeningAngle = 0.95f;
  m_contactParameters.m_distanceThreshold = 2.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int PhotonFragmentRemovalAlgorithm::GetMaxPasses() const { return m_nMaxPasses; }

//------------------------------------------------------------------------------------------------------------------------------------------

bool PhotonFragmentRemovalAlgorithm::IsCandidateDaughter(const Cluster* const pDaughterCluster) const {
  if (!pDaughterCluster->GetAssociatedTrackList().empty())
    return false;

  return (!m_useOnlyPhotonLikeDaughters || this->IsPhotonLike(pDaughterCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PhotonFragmentRemovalAlgorithm::IsCandidateParent(const Cluster* const pDaughterCluster,
                                                       const Cluster* const pParentCluster) const {
  if (!pParentCluster->GetAssociatedTrackList().empty())
    return false;

  if (pParentCluster->GetInnerPseudoLayer() > pDaughterCluster->GetInnerPseudoLayer() + m_innerLayerTolerance)
    return false;

  if (pDaughterCluster->GetInitialDirection().GetCosOpeningAngle(pParentCluster->GetInitialDirection()) <
      m_minCosOpeningAngle)
    return false;

  return pParentCluster->PassPhotonId(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PhotonFragmentRemovalAlgorithm::PostMergeAction(const Cluster* const pBestParentCluster) {
  PandoraContentApi::Cluster::Metadata metadata;
  metadata.m_particleId = PHOTON;

  return PandoraContentApi::Cluster::AlterMetadata(*this, pBestParentCluster, metadata);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PhotonFragmentRemovalAlgorithm::IsPhotonLike(const Cluster* const pDaughterCluster) const {
  if (pDaughterCluster->PassPhotonId(this->GetPandora()))
    return true;

  const ClusterFitResult& clusterFitResult(pDaughterCluster->GetFitToAllHitsResult());

  if ((PandoraContentApi::GetGeometry(*this)->GetHitTypeGranularity(pDaughterCluster->GetInnerLayerHitType()) <=
       FINE) &&
      (pDaughterCluster->GetInnerPseudoLayer() < m_photonLikeMaxInnerLayer) && (clusterFitResult.IsFitSuccessful()) &&
      (clusterFitResult.GetRadialDirectionCosine() > m_photonLikeMinDCosR) &&
      (pDaughterCluster->GetShowerProfileStart(this->GetPandora()) < m_photonLikeMaxShowerStart) &&
      (pDaughterCluster->GetShowerProfileDiscrepancy(this->GetPandora()) < m_photonLikeMaxProfileDiscrepancy)) {
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PhotonFragmentRemovalAlgorithm::PassesClusterContactCuts(const ClusterContact& clusterContact) const {
  if ((clusterContact.GetNContactLayers() > m_contactCutNLayers) ||
      (clusterContact.GetConeFraction1() > m_contactCutConeFraction1) ||
      (clusterContact.GetCloseHitFraction1() > m_contactCutCloseHitFraction1) ||
      (clusterContact.GetCloseHitFraction2() > m_contactCutCloseHitFraction2)) {
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PhotonFragmentRemovalAlgorithm::GetClusterMergingCandidates(const ContactMap& clusterContactMap,
                                                                       const Cluster*& pBestParentCluster,
                                                                       const Cluster*& pBestDaughterCluster) {
  float highestEvidence(m_minEvidence);
  float highestEvidenceParentEnergy(0.);

  ClusterList clusterList;
  for (const auto& mapEntry : clusterContactMap)
    clusterList.push_back(mapEntry.first);
  clusterList.sort(SortingHelper::SortClustersByNHits);

  for (const Cluster* const pDaughterCluster : clusterList) {
    const ContactVector& contactVector(clusterContactMap.at(pDaughterCluster));

    for (const ClusterContact& clusterContact : contactVector) {
      if (pDaughterCluster != clusterContact.GetDaughterCluster())
        throw StatusCodeException(STATUS_CODE_FAILURE);

      const float evidence(this->GetEvidenceForMerge(clusterContact));
      const float parentEnergy(clusterContact.GetParentCluster()->GetHadronicEnergy());

      if ((evidence > highestEvidence) ||
          ((evidence == highestEvidence) && (parentEnergy > highestEvidenceParentEnergy))) {
        highestEvidence = evidence;
        pBestDaughterCluster = pDaughterCluster;
        pBestParentCluster = clusterContact.GetParentCluster();
        highestEvidenceParentEnergy = parentEnergy;
      }
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float PhotonFragmentRemovalAlgorithm::GetEvidenceForMerge(const ClusterContact& clusterContact) const {
  // Calculate a measure of the evidence that the daughter candidate cluster is a fragment of the parent candidate
  // cluster:

  // 1. Layers in contact
  float contactEvidence(0.f);
  if ((clusterContact.GetNContactLayers() > m_contactEvidenceNLayers) &&
      (clusterContact.GetContactFraction() > m_contactEvidenceFraction)) {
    contactEvidence = clusterContact.GetContactFraction();
  }

  // 2. Cone extrapolation
  float coneEvidence(0.f);
  if (clusterContact.GetConeFraction1() > m_coneEvidenceFraction1) {
    coneEvidence = clusterContact.GetConeFraction1();
  }

  // 3. Distance of closest approach
  float distanceEvidence(0.f);
  if (clusterContact.GetDistanceToClosestHit() < m_distanceEvidence1) {
    distanceEvidence = (m_distanceEvidence1 - clusterContact.GetDistanceToClosestHit()) / m_distanceEvidence1d;
    distanceEvidence += m_distanceEvidenceCloseFraction1Multiplier * clusterContact.GetCloseHitFraction1();
    distanceEvidence += m_distanceEvidenceCloseFraction2Multiplier * clusterContact.GetCloseHitFraction2();
  }

  return ((m_contactWeight * contactEvidence) + (m_coneWeight * coneEvidence) + (m_distanceWeight * distanceEvidence));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PhotonFragmentRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadCommonSettings(xmlHandle));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "NMaxPasses", m_nMaxPasses));

  // Initial cluster candidate selection
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "InnerLayerTolerance", m_innerLayerTolerance));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinCosOpeningAngle", m_minCosOpeningAngle));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "UseOnlyPhotonLikeDaughters", m_useOnlyPhotonLikeDaughters));

  // Photon-like cuts
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "PhotonLikeMaxInnerLayer", m_photonLikeMaxInnerLayer));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "PhotonLikeMinDCosR", m_photonLikeMinDCosR));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "PhotonLikeMaxShowerStart", m_photonLikeMaxShowerStart));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "PhotonLikeMaxProfileDiscrepancy", m_photonLikeMaxProfileDiscrepancy));

  // Cluster contact cuts
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ContactCutNLayers", m_contactCutNLayers));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ContactCutConeFraction1", m_contactCutConeFraction1));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ContactCutCloseHitFraction1", m_contactCutCloseHitFraction1));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ContactCutCloseHitFraction2", m_contactCutCloseHitFraction2));

  // Total evidence: Contact evidence
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ContactEvidenceNLayers", m_contactEvidenceNLayers));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ContactEvidenceFraction", m_contactEvidenceFraction));

  // Cone evidence
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ConeEvidenceFraction1", m_coneEvidenceFraction1));

  // Distance of closest approach evidence
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "DistanceEvidence1", m_distanceEvidence1));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "DistanceEvidence1d", m_distanceEvidence1d));

  if (m_distanceEvidence1d < std::numeric_limits<float>::epsilon())
    return STATUS_CODE_INVALID_PARAMETER;

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "DistanceEvidenceCloseFraction1Multiplier",
                                                       m_distanceEvidenceCloseFraction1Multiplier));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "DistanceEvidenceCloseFraction2Multiplier",
                                                       m_distanceEvidenceCloseFraction2Multiplier));

  // Evidence weightings
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ContactWeight", m_contactWeight));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ConeWeight", m_coneWeight));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "DistanceWeight", m_distanceWeight));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinEvidence", m_minEvidence));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
