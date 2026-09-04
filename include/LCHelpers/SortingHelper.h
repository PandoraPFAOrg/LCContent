/**
 *  @file   LCContent/include/LCHelpers/SortingHelper.h
 *
 *  @brief  Header file for the sorting helper class.
 *
 *  $Log: $
 */
#ifndef LC_SORTING_HELPER_H
#define LC_SORTING_HELPER_H 1

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"

namespace lc_content {

/**
 *  @brief  SortingHelper class
 */
class SortingHelper {
public:
  /**
   *  @brief  Sort clusters by number of hits, associated track energy if no hits, then hadronic energy
   *
   *  @param  pLhs address of first cluster
   *  @param  pRhs address of second cluster
   */
  static bool SortClustersByNHits(const pandora::Cluster* const pLhs, const pandora::Cluster* const pRhs);

  /**
   *  @brief  Sort clusters by ascending inner layer, and by hadronic energy within a layer
   *
   *  @param  pLhs address of first cluster
   *  @param  pRhs address of second cluster
   */
  static bool SortClustersByInnerLayer(const pandora::Cluster* const pLhs, const pandora::Cluster* const pRhs);

  /**
   *  @brief  Sort pfos by descending energy
   *
   *  @param  pLhs address of first pfo
   *  @param  pRhs address of second pfo
   */
  static bool SortPfosByEnergy(const pandora::ParticleFlowObject* const pLhs,
                               const pandora::ParticleFlowObject* const pRhs);

  /**
   *  @brief  Sort tracks by ascending distance to cluster centroid at calorimeter
   *
   *  @param  pLhs address of first track
   *  @param  pRhs address of second track
   *  @param  pCluster address of cluster to calculate distance to
   */
  static bool SortTracksByDistance(const pandora::Track* a, const pandora::Track* b,
                                   const pandora::Cluster* const pCluster);
};

} // namespace lc_content

#endif // #ifndef LC_SORTING_HELPER_H
