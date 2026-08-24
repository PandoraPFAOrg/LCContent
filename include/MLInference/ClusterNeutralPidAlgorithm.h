/**
 *  @file   LCContent/include/MLInference/ClusterNeutralPidAlgorithm.h
 *
 *  @brief  Header file for the neutral-PID (sat/NH/photon) ONNX tagging algorithm.
 *
 *  Port of k4RecCalorimeter ClusterNeutralPidOnnx into a Pandora algorithm.  Runs a
 *  ParT multi-class model returning a softmax over {satellite, neutral-hadron, photon}
 *  for every cluster in the current list and writes the verdict as the cluster
 *  ParticleId:  photon -> pandora::PHOTON, NH -> pandora::K_LONG, satellite ->
 *  pandora::UNKNOWN_PARTICLE_TYPE.  Charged identification is left to the downstream
 *  track-cluster matching (the model was trained with charged primaries excluded).
 *
 *  Token construction matches photonId_hitcollect_multi.py / ClusterNeutralPidOnnx:
 *    own_ECAL  (<=BudgetOwnEcal)   cluster's own ECAL hits            type_id=0
 *    other_ECAL(<=BudgetOtherEcal) ECAL hits in cone, not in cluster  type_id=1
 *    HCAL      (<=BudgetHcal)      HCAL hits in cone                  type_id=2
 *    track     (<=BudgetTrack)     track impacts in cone             type_id=3
 *  Per-token features (5): (dTheta, dPhi, depth|FA, ln E|ln|p|, is_cher|FA), where
 *  FA = FEATURE_ABSENT (-1) marks a feature that does not apply to that token type and
 *  must stay -1 to match the trained model.
 */
#ifndef LC_CLUSTER_NEUTRAL_PID_ALGORITHM_H
#define LC_CLUSTER_NEUTRAL_PID_ALGORITHM_H 1

#include "Objects/CaloHit.h"
#include "Pandora/Algorithm.h"

#include "MLInference/OnnxSession.h"

#include <memory>
#include <string>
#include <vector>

namespace lc_content {

/**
 *  @brief  ClusterNeutralPidAlgorithm class
 */
class ClusterNeutralPidAlgorithm : public pandora::Algorithm {
public:
  /**
   *  @brief  Factory class for instantiating algorithm
   */
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm* CreateAlgorithm() const;
  };

  ClusterNeutralPidAlgorithm();
  ~ClusterNeutralPidAlgorithm();

private:
  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  // ---- token feature constants (must match the trained model) ----
  static constexpr float FEATURE_ABSENT = -1.f; ///< sentinel: feature N/A for this token type
  static constexpr int N_FEAT = 5;
  static constexpr int N_CLASSES = 3;
  static constexpr int TID_OWN = 0;
  static constexpr int TID_OTHER = 1;
  static constexpr int TID_HCAL = 2;
  static constexpr int TID_TRACK = 3;

  /**
   *  @brief  Pre-decoded calorimeter hit (Pandora-side)
   */
  struct Hit {
    float ux, uy, uz; ///< unit direction
    float th, ph;     ///< polar / azimuth
    float e;          ///< input (digi-channel) energy = CaloHit::GetInputEnergy()
    float depthFeat;  ///< ECAL: CaloHit::GetLayer(); HCAL token slot uses FEATURE_ABSENT
    float cherFeat;   ///< 1 if HitType==DRC_CHEREN, else 0
    bool isEcal;      ///< electromagnetic (ECAL) hit (GetHadronicEnergy() <= 0)
    const pandora::CaloHit* pCaloHit;
  };

  struct TrackImpact {
    float ux, uy, uz;
    float th, ph;
    float p;
  };

  pandora::StatusCode LoadModel();
  int Decide(float s0, float s1, float s2) const; ///< raw model scores [sat, NH, photon]

  // ---- configuration (XML) ----
  std::string m_modelPath; ///< path to the ParT multi-class ONNX file
  float m_w0;              ///< W0 for the log-weighted cluster barycenter (matches k4)
  float m_coneLateralMm;   ///< lateral cone [mm]: half-angle = atan(d / r_cluster)
  int m_budgetOwnEcal;     ///< max own-ECAL tokens per cluster
  int m_budgetOtherEcal;   ///< max other-ECAL tokens per cluster
  int m_budgetHcal;        ///< max HCAL tokens per cluster
  int m_budgetTrack;       ///< max track tokens per cluster
  float m_photonThreshold; ///< if <0: argmax(3); else photon if sm[photon] > threshold
  int m_maxTokens;         ///< sum of budgets

  // ---- ONNX session ----
  std::unique_ptr<OnnxSession> m_session; ///< ONNX session; null / invalid until the model is loaded
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm* ClusterNeutralPidAlgorithm::Factory::CreateAlgorithm() const {
  return new ClusterNeutralPidAlgorithm();
}

} // namespace lc_content

#endif // #ifndef LC_CLUSTER_NEUTRAL_PID_ALGORITHM_H
