/**
 *  @file   LCContent/include/MLInference/SatelliteAssignmentOnnxAlgorithm.h
 *
 *  @brief  Cluster-grouping satellite assignment via an ONNX pairwise-affinity model (SPLIT policy).
 *
 *  For each trackless, non-photon cluster (a satellite / neutral candidate) the ONNX model scores,
 *  for every in-cone neighbour cluster, the affinity P(same particle as the candidate).  SPLIT
 *  policy: the candidate is merged into the TRACKED cluster with the highest affinity if that
 *  affinity exceeds AffinityThreshold; otherwise it is kept as a neutral PFO.  There is NO
 *  neutral-neutral grouping and NO transitive union -- a candidate's charged/neutral fate depends
 *  only on its own DIRECT affinity to a tracked cluster (cascade-free).  Two tracked clusters are
 *  never merged (each track = one charged particle).  Isolated-hit merging is a downstream step.
 *
 *  The ONNX graph takes the candidate's token set -- own hits + in-cone neighbour hits + track
 *  tokens, tagged by type {own,other,track} and a pfo-index {0=candidate, 1..K=neighbours} -- and
 *  returns one affinity logit per neighbour pfo.  The pfo dimension is dynamic, so the number of
 *  neighbour clusters is never capped.
 */
#ifndef LC_SATELLITE_ASSIGNMENT_ONNX_ALGORITHM_H
#define LC_SATELLITE_ASSIGNMENT_ONNX_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/Cluster.h"
#include "Objects/CartesianVector.h"

#include "MLInference/OnnxSession.h"

#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace lc_content {

class SatelliteAssignmentOnnxAlgorithm : public pandora::Algorithm {
public:
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm *CreateAlgorithm() const override { return new SatelliteAssignmentOnnxAlgorithm(); }
  };

  SatelliteAssignmentOnnxAlgorithm() = default;
  ~SatelliteAssignmentOnnxAlgorithm() override;

private:
  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  // ---- token layout (must match the model's training tokenisation) ----
  static constexpr int N_FEAT  = 6;          ///< [dTheta, dPhi, depth, ln(E|p), isCher, isEcal]
  enum TokenType { TYPE_OWN = 0, TYPE_OTHER = 1, TYPE_TRACK = 2 };
  static constexpr float SENT = -1.0f;       ///< sentinel for absent depth/cher/ecal (track tokens)

  /// One calorimeter hit's per-token features (position reduced to angles).
  struct Hit {
    float theta = 0.f, phi = 0.f, energy = 0.f, depth = SENT;
    bool  isEcal = false, isCher = false;
  };

  /// Per-cluster cache: hits, energy-weighted centroid direction, and (if any) the track state.
  struct ClusterInfo {
    const pandora::Cluster  *pCluster = nullptr;
    pandora::CartesianVector centroidDir{0.f, 0.f, 0.f};
    std::vector<Hit>         hits;
    bool  hasTrack = false;                  ///< has a track reaching the calorimeter (any p) -> charged primary
    float trkTheta = 0.f, trkPhi = 0.f, trkP = 0.f;   ///< track state AT CALORIMETER (no projection)
  };

  /// One in-cone neighbour's model affinity to the candidate, tagged tracked/neutral.
  struct NeighbourAff {
    const pandora::Cluster *pCluster = nullptr;
    bool  hasTrack = false;                   ///< neighbour is a charged primary (has a calo-reaching track)
    float affinity = 0.f;                     ///< sigmoid(logit), P(same particle as candidate), in [0,1]
  };

  void LoadModel();
  static Hit ExtractHit(const pandora::CaloHit *pHit);

  /// Fill hits + energy-weighted centroid + best track state at the calo. False if no usable hits.
  bool BuildInfo(const pandora::Cluster *pCluster, ClusterInfo &info) const;

  /// Snapshot the current cluster list into per-cluster caches (hits, energy-weighted centroid, best
  /// track state).  Clusters with no usable hits are skipped.
  pandora::StatusCode BuildClusterCache(std::vector<ClusterInfo> &clusters) const;

  /// Select the merge candidates = available, trackless, non-photon clusters (as indices into
  /// `clusters`).  Tracked clusters are charged primaries (targets, never candidates); photons are
  /// left to the upstream photon logic.
  void CollectCandidates(const std::vector<ClusterInfo> &clusters, std::vector<int> &candidates) const;

  /// Build the candidate's token set (own + in-cone neighbour hits + track tokens), run the ONNX model
  /// ONCE, and return one affinity per in-cone neighbour cluster.  False on failure / no neighbour.
  bool ComputeNeighbourAffinities(const std::vector<ClusterInfo> &clusters, int candIndex,
                                  std::vector<NeighbourAff> &out) const;

  /// Pass 1 (frozen == deployed): plan each candidate's merge into its highest-affinity TRACKED
  /// neighbour if that affinity exceeds AffinityThreshold.  Fills chargedMerges (child -> parent).
  void PlanChargedMerges(const std::vector<ClusterInfo> &clusters, const std::vector<int> &candidates,
                         const std::map<int, std::vector<NeighbourAff>> &affinities,
                         std::map<const pandora::Cluster *, const pandora::Cluster *> &chargedMerges) const;

  /// Pass 2 (opt-in): union-find over neutral<->neutral affinities exceeding NeutralAffinityThreshold,
  /// among the clusters that stayed neutral (charged assignment frozen).  Fills neutralGroups; each
  /// group's representative (the biggest cluster, merged INTO) is placed first.
  void PlanNeutralGroups(const std::vector<ClusterInfo> &clusters, const std::vector<int> &candidates,
                         const std::map<int, std::vector<NeighbourAff>> &affinities,
                         const std::set<const pandora::Cluster *> &becameCharged,
                         std::vector<std::vector<const pandora::Cluster *>> &neutralGroups) const;

  /// Execute the planned charged merges (candidate cluster -> its tracked parent).
  pandora::StatusCode ApplyChargedMerges(
      const std::map<const pandora::Cluster *, const pandora::Cluster *> &merges) const;

  /// Execute the planned neutral groups (merge every other member of each component into its first =
  /// representative cluster).
  pandora::StatusCode ApplyNeutralGroups(
      const std::vector<std::vector<const pandora::Cluster *>> &groups) const;

  std::string m_modelPath;                   ///< ONNX model path (required)
  float m_affinityThreshold = 0.50f;         ///< merge to charged if best tracked affinity > this
  bool  m_groupNeutral = true;               ///< 2nd pass: union-find neutral<->neutral grouping (on by default)
  float m_neutralAffinityThreshold = 0.72f;  ///< union two neutral clusters if their affinity > this
  float m_coneHalfAngle     = 0.4f;          ///< rad, token/neighbour cone half-angle
  float m_coneCos           = std::cos(0.4f);
  int   m_budgetOwn         = 512;           ///< top-K own hits by energy
  int   m_budgetNbr         = 64;            ///< top-K hits per neighbour cluster by energy (no cap on #clusters)

  std::unique_ptr<OnnxSession> m_session;    ///< ONNX session; null / invalid until the model is loaded
};

} // namespace lc_content

#endif // LC_SATELLITE_ASSIGNMENT_ONNX_ALGORITHM_H
