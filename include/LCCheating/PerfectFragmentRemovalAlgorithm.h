/**
 *  @file   LCContent/include/LCCheating/PerfectFragmentRemovalAlgorithm.h
 *
 *  @brief  Header file for the perfect fragment removal algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_PERFECT_FRAGMENT_REMOVAL_ALGORITHM_H
#define LC_PERFECT_FRAGMENT_REMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief PerfectFragmentRemovalAlgorithm class
 */
class PerfectFragmentRemovalAlgorithm : public pandora::Algorithm {
public:
  /**
   *  @brief Default constructor
   */
  PerfectFragmentRemovalAlgorithm();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  typedef std::map<const pandora::MCParticle*, const pandora::Cluster*> MCParticleToClusterMap;

  bool m_shouldMergeChargedClusters; ///< Whether to merge charged clusters sharing same mc particle (otherwise use only
                                     ///< highest E charged cluster)
};

} // namespace lc_content

#endif // #ifndef LC_PERFECT_FRAGMENT_REMOVAL_ALGORITHM_H
