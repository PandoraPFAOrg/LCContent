#ifndef IdeaPfoCreationAlgorithm_h
#define IdeaPfoCreationAlgorithm_h 1

#include "Pandora/Algorithm.h"
#include "Api/PandoraContentApi.h"
#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"
#include "Objects/Track.h"

namespace lc_content {

class IdeaPfoCreationAlgorithm : public pandora::Algorithm {
public:
  IdeaPfoCreationAlgorithm();
  ~IdeaPfoCreationAlgorithm()=default;

  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm* CreateAlgorithm() const;
  };

private:
  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  pandora::StatusCode CreateElectronCandidates(const pandora::ClusterList& aList) const;
  pandora::StatusCode CreateChargedHadronCandidates(const pandora::ClusterList& aList) const;
  pandora::StatusCode CreatePhotonCandidates(const pandora::ClusterList& aList) const;
  pandora::StatusCode CreateNeutralHadronCandidates(const pandora::ClusterList& aList) const;
  pandora::StatusCode CreatePfoFromTrack(const pandora::TrackList* aList) const;

  pandora::StatusCode GetDualReadoutEnergy(const pandora::Cluster* aClus, float& energy) const;
  pandora::StatusCode IsEmShower(const pandora::Cluster* aClus, bool& isEmShower) const;
  float EstimateScintEnergy(const pandora::Cluster* aClus) const;

  const pandora::Track* FindBestAssociatedTrack(const pandora::Cluster* aClus) const;

  /// Calorimeter energy resolution: sigma = stoch/sqrt(p) + const (floored at m_sigmaFloor)
  float CaloSigma(float p, bool isEm) const;

  std::string m_outputPfoListName;

  float m_nSigma         = 2.0f;   ///< Threshold in sigma units for the E_DR/p gate (validated optimum n=2)
  float m_stochasticHad  = 0.30f;  ///< Hadronic stochastic term
  float m_constantHad    = 0.01f;  ///< Hadronic constant term
  float m_stochasticEm   = 0.02f;  ///< EM stochastic term
  float m_constantEm     = 0.005f; ///< EM constant term
  float m_ptCut          = 0.6f;

  // ---- regional-excess significance gate for the charged-hadron energy ----
  float m_coneOpeningAngle = 0.4f;  ///< Regional cone opening angle [rad] for the E/p significance
  float m_significanceK    = 2.0f;  ///< gate: emit iff  S = (sumEdr-sumP)/(sigma_c(sumEdr)*sumEdr) > k  (validated k=2)

  // non-configurable safe guard
  float m_sigmaFloor     = 0.001f;  ///< Minimum sigma value

  /// Flat hadronic-response calibration multiplying the E_DR neutral-hadron PFO energy
  /// (default 1.0 = no scaling).  XML "NeutralHadEnergyScale".
  float m_neutralHadScale = 1.0f;
}; // class IdeaPfoCreationAlgorithm

inline pandora::Algorithm* IdeaPfoCreationAlgorithm::Factory::CreateAlgorithm() const {
  return new IdeaPfoCreationAlgorithm();
}

} // namespace lc_content

#endif // #ifndef IdeaPfoCreationAlgorithm_h
