/**
 *  @file   LCContent/include/LCTopologicalAssociation/BIBHitMergingAlgorithm.h
 *
 *  @brief  Header file for the BIB hit merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_BIB_HIT_MERGING_ALGORITHM_H
#define LC_BIB_HIT_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <utility>

namespace lc_content {

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

/**
 *  @brief  BIBHitMergingAlgorithm class
 *
 *  Controlled recovery of calo hits flagged isPossibleBIB: a flagged, available hit is
 *  re-attached (as an isolated hit) to a cluster only when it is locally continuous with
 *  that cluster - at least MinClusterNeighbours of the hit's neighbours (same transverse
 *  cylinder as the flag definition) are non-isolated member hits of the cluster - and the
 *  cluster is an eligible host (track-associated, or above MinHostClusterEnergy), the hit
 *  lies within the cluster's pseudolayer span (+- LayerSpanTolerance), and the nearest
 *  member hit is within MaxRecombinationDistance. Shower halo touches its shower and
 *  passes; BIB carpet hits near a cluster but not woven into it do not. Attached hits do
 *  not become members, so recovery cannot chain outward across the BIB carpet.
 *
 *  Optional timing gate (TimingCutEnabled, default off): additionally require the hit
 *  time to be within MaxTimeDifference of the host cluster's energy-weighted mean hit
 *  time. Intended to run after track-cluster association, so track matches are known.
 */
class BIBHitMergingAlgorithm : public pandora::Algorithm {
public:
  typedef KDTreeLinkerAlgo<std::pair<const pandora::CaloHit*, unsigned int>, 3> HitKDTree3D;
  typedef KDTreeNodeInfoT<std::pair<const pandora::CaloHit*, unsigned int>, 3> HitKDNode3D;

  /**
   *  @brief Default constructor
   */
  BIBHitMergingAlgorithm();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  bool m_shouldUseCurrentClusterList; ///< Whether to use clusters from the current list in the algorithm
  pandora::StringVector m_additionalClusterListNames; ///< Additional cluster lists from which to consider clusters

  float m_minHostClusterEnergy; ///< Host eligibility: min cluster hadronic energy for clusters without an associated
                                ///< track, units GeV
  float m_maxRecombinationDistance; ///< Max distance between the hit and the nearest non-isolated member hit of the
                                    ///< host, units mm
  unsigned int
      m_minClusterNeighbours; ///< Min number of the hit's cylinder neighbours that must belong to the host cluster
  unsigned int
      m_neighbourNLayers; ///< Pseudolayer window of the neighbour cylinder (same convention as the flag definition)
  float m_neighbourCutDistanceSquared;   ///< Transverse distance of the neighbour cylinder, stored squared, units mm^2
                                         ///< (XML value NeighbourCutDistance is in mm)
  float m_neighbourMaxSeparationSquared; ///< Max 3D separation of cylinder neighbours, stored squared, units mm^2 (XML
                                         ///< value NeighbourMaxSeparation is in mm)
  float m_searchSafetyFactor; ///< Safety factor applied to the neighbour cut distance to define the kd-tree search box
  unsigned int m_layerSpanTolerance; ///< Allowed pseudolayer excursion of the hit beyond the host cluster's [inner,
                                     ///< outer] layer span

  bool m_timingCutEnabled;   ///< Whether to require time compatibility between the hit and the host cluster
  float m_maxTimeDifference; ///< Max |hit time - host cluster energy-weighted mean hit time|, units ns
};

} // namespace lc_content

#endif // #ifndef LC_BIB_HIT_MERGING_ALGORITHM_H
