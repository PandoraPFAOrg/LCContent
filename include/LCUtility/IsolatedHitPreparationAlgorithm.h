#ifndef IsolatedHitPreparationAlgorithm_h
#define IsolatedHitPreparationAlgorithm_h 1

#include "Pandora/Algorithm.h"

namespace lc_content {

class IsolatedHitPreparationAlgorithm : public pandora::Algorithm {
public:
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm* CreateAlgorithm() const override;
  };

  IsolatedHitPreparationAlgorithm() = default;

private:
  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  bool m_excludeEcal = false;         ///< whether to exclude ECAL hits from isolated hit flagging (default: false)
  bool m_excludeHcal = false;         ///< whether to exclude HCAL hits from isolated hit flagging (default: false)
  float m_scintEnergyThreshold = 0.f; ///< energy cut (GeV) applied to scint hits (and cher if cher cut = 0)
  float m_cherEnergyThreshold = 0.f;  ///< if > 0, energy cut (GeV) applied to cherenkov hits instead
};

} // namespace lc_content

#endif // IsolatedHitPreparationAlgorithm_h
