/**
 *  @file   LCContent/include/LCUtility/ClusterPreparationAlgorithm.h
 *
 *  @brief  Header file for the cluster preparation algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CLUSTER_PREPARATION_ALGORITHM_H
#define LC_CLUSTER_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief  ClusterPreparationAlgorithm class
 */
class ClusterPreparationAlgorithm : public pandora::Algorithm {
private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  pandora::StringVector m_candidateListNames; ///< The list of cluster list names to use
  std::string m_mergedCandidateListName;      ///< The name of the merged candidate list name
};

} // namespace lc_content

#endif // #ifndef LC_CLUSTER_PREPARATION_ALGORITHM_H
