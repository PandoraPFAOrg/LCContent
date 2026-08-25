/**
 *  @file   LCContent/include/LCFragmentRemoval/PhotonFragmentRemovalAlgorithm.h
 *
 *  @brief  Header file for the photon fragment removal algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_PHOTON_FRAGMENT_REMOVAL_ALGORITHM_H
#define LC_PHOTON_FRAGMENT_REMOVAL_ALGORITHM_H 1

#include "LCFragmentRemoval/FragmentRemovalBaseAlgorithm.h"

#include "LCHelpers/FragmentRemovalHelper.h"

namespace lc_content {

/**
 *  @brief  PhotonFragmentRemovalAlgorithm class
 */
class PhotonFragmentRemovalAlgorithm : public FragmentRemovalBaseAlgorithm<ClusterContact> {
public:
  /**
   *  @brief Default constructor
   */
  PhotonFragmentRemovalAlgorithm();

private:
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  // Hooks supplying what is specific to this algorithm; see FragmentRemovalBaseAlgorithm for what each is asked
  unsigned int GetMaxPasses() const;
  bool IsCandidateDaughter(const pandora::Cluster* const pDaughterCluster) const;
  bool IsCandidateParent(const pandora::Cluster* const pDaughterCluster,
                         const pandora::Cluster* const pParentCluster) const;
  bool PassesClusterContactCuts(const ClusterContact& clusterContact) const;

  pandora::StatusCode GetClusterMergingCandidates(const ContactMap& clusterContactMap,
                                                  const pandora::Cluster*& pBestParentCluster,
                                                  const pandora::Cluster*& pBestDaughterCluster);

  /**
   *  @brief  Relabel the enlarged parent cluster as a photon
   *
   *  @param  pBestParentCluster address of the parent cluster of the merge just applied
   */
  pandora::StatusCode PostMergeAction(const pandora::Cluster* const pBestParentCluster);

  /**
   *  @brief  Whether candidate daughter cluster can be considered as photon-like
   *
   *  @param  pDaughterCluster address of the candidate daughter cluster
   *
   *  @return boolean
   */
  bool IsPhotonLike(const pandora::Cluster* const pDaughterCluster) const;

  /**
   *  @brief  Get a measure of the evidence for merging the parent and daughter candidate clusters
   *
   *  @param  clusterContact the cluster contact details for parent/daughter candidate merge
   *
   *  @return the evidence
   */
  float GetEvidenceForMerge(const ClusterContact& clusterContact) const;

  unsigned int m_nMaxPasses; ///< Maximum number of passes over cluster contact information

  unsigned int m_innerLayerTolerance; ///< Max number of layers by which daughter can exceed parent inner layer
  float m_minCosOpeningAngle;         ///< Min cos opening angle between candidate cluster initial directions

  bool m_useOnlyPhotonLikeDaughters; ///< Whether to skip photon-like checks for daughter cluster

  unsigned int m_photonLikeMaxInnerLayer;  ///< Max inner layer to identify daughter cluster as photon-like
  float m_photonLikeMinDCosR;              ///< Max radial direction cosine to identify daughter as photon-like
  float m_photonLikeMaxShowerStart;        ///< Max shower profile start to identify daughter as photon-like
  float m_photonLikeMaxProfileDiscrepancy; ///< Max shower profile discrepancy to identify daughter as photon-like

  unsigned int m_contactCutNLayers;    ///< Number of contact layers to store cluster contact info
  float m_contactCutConeFraction1;     ///< Cone fraction 1 value to store cluster contact info
  float m_contactCutCloseHitFraction1; ///< Close hit fraction 1 value to store cluster contact info
  float m_contactCutCloseHitFraction2; ///< Close hit fraction 2 value to store cluster contact info

  unsigned int m_contactEvidenceNLayers;            ///< Contact layers required for contact evidence contribution
  float m_contactEvidenceFraction;                  ///< Contact fraction required for contact evidence contribution
  float m_coneEvidenceFraction1;                    ///< Cone fraction 1 value required for cone evidence contribution
  float m_distanceEvidence1;                        ///< Offset for distance evidence contribution 1
  float m_distanceEvidence1d;                       ///< Denominator for distance evidence contribution 1
  float m_distanceEvidenceCloseFraction1Multiplier; ///< Distance evidence multiplier for close hit fraction 1
  float m_distanceEvidenceCloseFraction2Multiplier; ///< Distance evidence multiplier for close hit fraction 2

  float m_contactWeight;  ///< Weight for layers in contact evidence
  float m_coneWeight;     ///< Weight for cone extrapolation evidence
  float m_distanceWeight; ///< Weight for distance of closest approach evidence

  float m_minEvidence; ///< Min evidence before parent/daughter candidates can be merged
};

} // namespace lc_content

#endif // #ifndef LC_PHOTON_FRAGMENT_REMOVAL_ALGORITHM_H
