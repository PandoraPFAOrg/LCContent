/**
 *  @file   LCContent/include/LCFragmentRemoval/NeutralFragmentRemovalAlgorithm.h
 *
 *  @brief  Header file for the neutral fragment removal algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_NEUTRAL_FRAGMENT_REMOVAL_ALGORITHM_H
#define LC_NEUTRAL_FRAGMENT_REMOVAL_ALGORITHM_H 1

#include "LCFragmentRemoval/FragmentRemovalBaseAlgorithm.h"

#include "LCHelpers/ClusterProximityHelper.h"
#include "LCHelpers/FragmentRemovalHelper.h"

namespace lc_content {

/**
 *  @brief  NeutralClusterContact class, describing the interactions and proximity between parent and daughter candidate
 * clusters
 */
class NeutralClusterContact : public ClusterContact {
public:
  /**
   *  @brief  Parameters class
   */
  class Parameters : public ClusterContact::Parameters {
  public:
    float m_coneCosineHalfAngle2; ///< Cosine half angle for second cone comparison in cluster contact object
    float m_coneCosineHalfAngle3; ///< Cosine half angle for third cone comparison in cluster contact object
  };

  /**
   *  @brief  Constructor
   *
   *  @param  pandora the associated pandora instance
   *  @param  pDaughterCluster address of the daughter candidate cluster
   *  @param  pParentCluster address of the parent candidate cluster
   *  @param  parameters the cluster contact parameters
   *  @param  contactCache cache of the per-cluster bounding boxes the hit comparison uses
   */
  NeutralClusterContact(const pandora::Pandora& pandora, const pandora::Cluster* const pDaughterCluster,
                        const pandora::Cluster* const pParentCluster, const Parameters& parameters,
                        ClusterContactCache& contactCache);

  /**
   *  @brief  Get the fraction of daughter hits that lie within specified cone 2 along parent direction
   *
   *  @return The daughter cone fraction
   */
  float GetConeFraction2() const;

  /**
   *  @brief  Get the fraction of daughter hits that lie within specified cone 3 along parent direction
   *
   *  @return The daughter cone fraction
   */
  float GetConeFraction3() const;

private:
  float m_coneFraction2; ///< Fraction of daughter hits that lie within specified cone 2 along parent direction
  float m_coneFraction3; ///< Fraction of daughter hits that lie within specified cone 3 along parent direction
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  NeutralFragmentRemovalAlgorithm class
 */
class NeutralFragmentRemovalAlgorithm : public FragmentRemovalBaseAlgorithm<NeutralClusterContact> {
public:
  /**
   *  @brief Default constructor
   */
  NeutralFragmentRemovalAlgorithm();

private:
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  // Hooks supplying what is specific to this algorithm; see FragmentRemovalBaseAlgorithm for what each is asked
  unsigned int GetMaxPasses() const;
  bool IsCandidateDaughter(const pandora::Cluster* const pDaughterCluster) const;
  bool IsCandidateParent(const pandora::Cluster* const pDaughterCluster,
                         const pandora::Cluster* const pParentCluster) const;
  bool PassesClusterContactCuts(const NeutralClusterContact& neutralClusterContact) const;

  pandora::StatusCode GetClusterMergingCandidates(const ContactMap& neutralClusterContactMap,
                                                  const pandora::Cluster*& pBestParentCluster,
                                                  const pandora::Cluster*& pBestDaughterCluster);

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
   *  @param  neutralClusterContact the cluster contact details for parent/daughter candidate merge
   *
   *  @return the evidence
   */
  float GetEvidenceForMerge(const NeutralClusterContact& neutralClusterContact) const;

  unsigned int m_nMaxPasses; ///< Maximum number of passes over cluster contact information

  unsigned int m_photonLikeMaxInnerLayer;  ///< Max inner layer to identify daughter cluster as photon-like
  float m_photonLikeMinDCosR;              ///< Max radial direction cosine to identify daughter as photon-like
  float m_photonLikeMaxShowerStart;        ///< Max shower profile start to identify daughter as photon-like
  float m_photonLikeMaxProfileDiscrepancy; ///< Max shower profile discrepancy to identify daughter as photon-like

  unsigned int m_contactCutNLayers;          ///< Number of contact layers to store cluster contact info
  float m_contactCutConeFraction1;           ///< Cone fraction 1 value to store cluster contact info
  float m_contactCutCloseHitFraction1;       ///< Close hit fraction 1 value to store cluster contact info
  float m_contactCutCloseHitFraction2;       ///< Close hit fraction 2 value to store cluster contact info
  float m_contactCutNearbyDistance;          ///< Distance between closest hits to mark clusters as nearby
  float m_contactCutNearbyCloseHitFraction2; ///< Close hit fraction 2 in nearby hits to store cluster contact info

  unsigned int m_contactEvidenceNLayers1; ///< Contact evidence n layers cut 1
  unsigned int m_contactEvidenceNLayers2; ///< Contact evidence n layers cut 2
  unsigned int m_contactEvidenceNLayers3; ///< Contact evidence n layers cut 3
  float m_contactEvidence1;               ///< Contact evidence contribution 1
  float m_contactEvidence2;               ///< Contact evidence contribution 2
  float m_contactEvidence3;               ///< Contact evidence contribution 3

  float m_coneEvidenceFraction1;                 ///< Cone fraction 1 value required for cone evidence contribution
  float m_coneEvidenceFineGranularityMultiplier; ///< Cone evidence multiplier for fine granularity daughter clusters

  float m_distanceEvidence1;                        ///< Offset for distance evidence contribution 1
  float m_distanceEvidence1d;                       ///< Denominator for distance evidence contribution 1
  float m_distanceEvidenceCloseFraction1Multiplier; ///< Distance evidence multiplier for close hit fraction 1
  float m_distanceEvidenceCloseFraction2Multiplier; ///< Distance evidence multiplier for close hit fraction 2

  float m_contactWeight;  ///< Weight for layers in contact evidence
  float m_coneWeight;     ///< Weight for cone extrapolation evidence
  float m_distanceWeight; ///< Weight for distance of closest approach evidence

  float m_minEvidence; ///< Min evidence before parent/daughter candidates can be merged
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NeutralClusterContact::GetConeFraction2() const { return m_coneFraction2; }

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NeutralClusterContact::GetConeFraction3() const { return m_coneFraction3; }

} // namespace lc_content

#endif // #ifndef LC_NEUTRAL_FRAGMENT_REMOVAL_ALGORITHM_H
