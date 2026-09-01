/**
 *  @file   LCContent/include/LCMonitoring/EfficiencyMonitoringAlgorithm.h
 *
 *  @brief  Header file for the efficiency monitoring algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_EFFICIENCY_MONITORING_ALGORITHM_H
#define LC_EFFICIENCY_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief EfficiencyMonitoringAlgorithm class
 */
class EfficiencyMonitoringAlgorithm : public pandora::Algorithm {
public:
  /**
   *  @brief Default constructor
   */
  EfficiencyMonitoringAlgorithm();

  /**
   *  @brief  Destructor
   */
  ~EfficiencyMonitoringAlgorithm();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  std::string m_monitoringFileName; ///< The name of the file in which to save the monitoring tree
  float m_mcThresholdEnergy;        ///< MC particle threshold energy, units GeV
};

} // namespace lc_content

#endif // #ifndef LC_EFFICIENCY_MONITORING_ALGORITHM_H
