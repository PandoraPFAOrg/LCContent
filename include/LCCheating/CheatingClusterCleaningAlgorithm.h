/**
 *  @file   LCContent/include/LCCheating/CheatingClusterCleaningAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CHEATING_CLUSTER_CLEANING_ALGORITHM_H
#define LC_CHEATING_CLUSTER_CLEANING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Pandora/PandoraInternal.h"

namespace lc_content {

/**
 *  @brief CheatingClusterCleaningAlgorithm class
 */
class CheatingClusterCleaningAlgorithm : public pandora::Algorithm {
private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lc_content

#endif // #ifndef LC_CHEATING_CLUSTER_CLEANING_ALGORITHM_H
