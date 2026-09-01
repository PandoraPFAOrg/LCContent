/**
 *  @file   LCContent/include/LCUtility/TrainingSoftwareCompensation.h
 *
 *  @brief  Header file for the training software compensation algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_TRAINING_SOFTWARE_COMPENSATION_H
#define LC_TRAINING_SOFTWARE_COMPENSATION_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief  TrainingSoftwareCompensation class
 */
class TrainingSoftwareCompensation : public pandora::Algorithm {
public:
  /**
   *  @brief  Default constructor
   */
  TrainingSoftwareCompensation();

  /**
   *  @breif Destructor
   */
  ~TrainingSoftwareCompensation();

private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  std::string m_myRootFileName; ///< Output root file for training of software compensation weights
  std::string
      m_trainingTreeName; ///< Name of the TTree in root file produced for training of software compensation weights
};

} // namespace lc_content

#endif // #ifndef LC_TRAINING_SOFTWARE_COMPENSATION_H
