/**
 *  @file   LCContent/include/LCUtility/PfoPreparationAlgorithm.h
 *
 *  @brief  Header file for the pfo preparation algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_PFO_PREPARATION_ALGORITHM_H
#define LC_PFO_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief  PfoPreparationAlgorithm class
 */
class PfoPreparationAlgorithm : public pandora::Algorithm {
private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  pandora::StringVector m_candidateListNames; ///< The list of pfo list names to use
  std::string m_mergedCandidateListName;      ///< The name of the merged candidate list name
};

} // namespace lc_content

#endif // #ifndef LC_PFO_PREPARATION_ALGORITHM_H
