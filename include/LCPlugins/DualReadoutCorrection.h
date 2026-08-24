/**
 *  @file   LCContent/include/LCPlugins/DualReadoutCorrection.h
 * 
 *  @brief  Header file for the dual-readout correction plugin algorithm class.
 * 
 *  $Log: $
 */
#ifndef DUALREADOUT_CORRECTION_H 
#define DUALREADOUT_CORRECTION_H 1

#include "Plugins/EnergyCorrectionsPlugin.h"

namespace lc_content {
/**
  *   @brief  Simple struct to hold dual readout energy information
  */
struct DualReadoutEnergy {
  float ecalS = 0.f;
  float ecalC = 0.f;
  float hcalS = 0.f;
  float hcalC = 0.f;
  float energyEcalCorrected = 0.f;
  float energyHcalCorrected = 0.f;

  // for particle ID
  // Note: SCEP cal doesn't have Cherenkov ch. at depth 0
  std::map<int, float> ecalS_byDepth;
  std::map<int, float> ecalC_byDepth;
};

/**
  *   @brief  A helper function to apply dual-readout correction to a cluster
  *           and fill the DualReadoutEnergy struct
  */
DualReadoutEnergy RunDualReadoutCorrection(const pandora::Cluster* const aClus,
                                           const float chiEcal, const float chiHcal);
/**
 *  @brief  DualReadoutCorrection class. 
 */
class DualReadoutCorrection : public pandora::EnergyCorrectionPlugin {
public:
  /**
  *  @brief  Constructor with input parameters
  *
  *  @param  parameters the input parameters 
  */
  DualReadoutCorrection()=default;

  pandora::StatusCode MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const override;

  /**
  *  @brief  Access the dual-readout correction factors, so that consumers of the corrected energy
  *          -- e.g. the error propagation in the edm4hep pfo creator -- use the same chi that was
  *          applied here instead of keeping a second copy of the values.
  */
  float GetChiEcal() const { return m_chiEcal; }
  float GetChiHcal() const { return m_chiHcal; }

private:
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  float m_chiEcal = 0.41f; // ECAL dual-readout correction factor
  float m_chiHcal = 0.31f; // HCAL dual-readout correction factor

}; // class DualReadoutCorrection
} // namespace lc_content

#endif // #ifndef DUALREADOUT_CORRECTION_H