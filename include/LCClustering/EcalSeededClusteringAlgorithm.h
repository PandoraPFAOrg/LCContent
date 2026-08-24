/**
 *  @file   LCContent/include/LCClustering/EcalSeededClusteringAlgorithm.h
 *
 *  @brief  Header for the ECAL-seeded HCAL clustering algorithm.
 *
 *  Grows HCAL clusters from ECAL cluster seeds, with BFS topological
 *  connectivity in (theta, phi) and opening-angle contest resolution on
 *  contested hits.  Unseeded HCAL clusters (late-developing neutral
 *  hadrons with no ECAL precursor) are also formed in a second pass.
 */

#ifndef ECAL_SEEDED_CLUSTERING_ALGORITHM_H
#define ECAL_SEEDED_CLUSTERING_ALGORITHM_H 1

#include "Objects/CartesianVector.h"
#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <vector>

namespace lc_content {

class EcalSeededClusteringAlgorithm : public pandora::Algorithm {
public:
  // ------------------------------------------------------------------ //
  //  Factory
  // ------------------------------------------------------------------ //
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm* CreateAlgorithm() const override { return new EcalSeededClusteringAlgorithm(); }
  };

  EcalSeededClusteringAlgorithm() = default;

private:
  // ================================================================== //
  //  Internal data structures
  // ================================================================== //
  struct HitCache;
  typedef std::map<int, std::vector<const HitCache*>> SeedCircleMap;

  /**
   *  @brief  Compact per-hit cache populated once at startup.
   *          Storing theta/phi avoids recomputing them in the inner loops.
   */
  struct HitCache {
    const pandora::CaloHit* pHit = nullptr;
    pandora::CartesianVector pos{0.f, 0.f, 0.f};
    float energy = 0.f; ///< hadronic energy
    float theta = 0.f;
    float phi = 0.f;
    bool isCheren = false;
  };

  /**
   *  @brief  Common state shared by all cluster types:
   *          the assigned hit list, the DR-corrected energy, and the
   *          log-weighted barycenter.  Updated by ComputeClusterState().
   */
  struct ClusterState {
    std::vector<const HitCache*> hits;                  ///< all assigned hits (pool owns storage)
    std::vector<const HitCache*> vicinity;              ///< candidate hits for BFS growth
    float energy = 0.f;                                 ///< DR-corrected running energy
    pandora::CartesianVector barycenter{0.f, 0.f, 0.f}; ///< log-weighted barycenter
  };

  /**
   *  @brief  One entry per surviving ECAL cluster seed.
   */
  struct SeededCluster : public ClusterState {
    const pandora::Cluster* pEcalCluster = nullptr; ///< target for AddToCluster
    bool isCharged = false;
    const pandora::Track* pBestTrack = nullptr; ///< highest-p track if charged

    /// Impact point on the HCAL inner surface (mm).
    pandora::CartesianVector hcalImpactPoint{0.f, 0.f, 0.f};
    float thetaSeed = 0.f;
    float phiSeed = 0.f;

    bool alive = true;
  };

  /**
   *  @brief  One entry per unseeded HCAL cluster (Phase 5).
   *          All state is inherited from ClusterState.
   */
  struct UnseededCluster : public ClusterState {};

  /// Contest-resolution policy for the shared Grow() engine.
  enum class GrowStrategy { ResolveByDist, MergeClusters };

  // ================================================================== //
  //  Phase methods
  // ================================================================== //

  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  /**
   *  @brief  Phase 1a - collect ECAL cluster seeds, filling SeededCluster.
   *          Photon clusters are skipped.  Charged seeds below
   *          m_pThres are dropped.
   */
  pandora::StatusCode SelectSeeds(const pandora::ClusterList& ecalClusters, std::vector<SeededCluster>& seeds) const;

  /**
   *  @brief  Phase 1b - collect all available HCAL hits into a flat cache.
   *          ECAL hits are identified by hadronic energy == 0 and skipped.
   *          Hits below m_hardThreshold are also skipped (scint only).
   */
  void BuildHcalPool(const pandora::CaloHitList& allHits, std::vector<HitCache>& hcalPool) const;

  /**
   *  @brief  Phase 2a - find the (theta, phi) impact point on the HCAL
   *          inner surface for each seed.
   *          Charged seeds: helix extrapolation from track state at ECAL.
   *          Neutral seeds:  radial projection of ECAL centroid.
   */
  pandora::StatusCode ExtrapolateToHcalSurface(SeededCluster& seed) const;

  /**
   *  @brief  Phase 2b - for each seed, collect the subset of hcalPool
   *          hits within m_maxSearchDeltaR of the seed impact point.
   */
  void BuildVicinities(std::vector<SeededCluster>& seeds, const std::vector<HitCache>& hcalPool) const;

  /**
   *  @brief  Phase 3a - for each seed, collect all vicinity hits within
   *          m_seedCircleRadius of the seed impact point that pass
   *          m_growThreshold.  Seeds whose scint sum falls below
   *          m_rhoSeed are killed here.
   */
  void FormSeedCircles(std::vector<SeededCluster>& seeds, SeedCircleMap& circles) const;

  /**
   *  @brief  Phase 3b - resolve hits that appear in more than one seed
   *          circle using the opening-angle distance against each seed's
   *          impact point.  Losing seeds keep the remaining
   *          hits (no second rho_seed check).
   */
  void ResolveSeedCircleContests(std::vector<SeededCluster>& seeds, SeedCircleMap& circles) const;

  /**
   *  @brief  Shared BFS growth engine used by both Phase 4 (seeded) and
   *          Phase 5 (unseeded).
   *
   *  Each ClusterState carries its own vicinity; Grow() only searches within
   *  that vicinity, avoiding the cost of scanning the full hit pool.
   *  For MergeClusters, when two clusters are merged their vicinities are
   *  unioned so the merged cluster can grow into both original regions.
   *
   *  ResolveByDist - contested hits go to the cluster with the lowest
   *                  opening-angle distance; barycenter/energy updated each layer.
   *  MergeClusters - contested hits cause the two competing clusters to be
   *                  merged (union-find); the merged cluster accumulates both.
   *
   *  @param strategy     Contest-resolution policy.
   *  @param clusters     Initial clusters; each must have vicinity populated.
   *  @param assignedSet  Already-assigned hits (blocked from re-assignment); updated in place.
   *  @return             Grown clusters (folded for MergeClusters).
   */
  std::vector<ClusterState> Grow(GrowStrategy strategy, const std::vector<ClusterState>& clusters,
                                 pandora::CaloHitSet& assignedSet) const;

  /**
   *  @brief  Phases 2-4 combined - extrapolate seeds to HCAL surface,
   *          build vicinities, form and contest seed circles, then BFS-grow.
   *          On return, assignedSet contains all HCAL hits claimed by seeds.
   */
  pandora::StatusCode FormSeededClusters(std::vector<SeededCluster>& seeds, const std::vector<HitCache>& hcalPool,
                                         pandora::CaloHitSet& assignedSet) const;

  /**
   *  @brief  Phase 5 - cluster leftover HCAL hits into unseeded clusters.
   *          Candidate seeds are hits whose local scint circle energy
   *          exceeds m_rhoUnseeded.  BFS growth with merge-on-contest
   *          via union-find.
   */
  pandora::StatusCode FormUnseededClusters(std::vector<UnseededCluster>& unseeded,
                                           const std::vector<HitCache>& hcalPool,
                                           pandora::CaloHitSet& assignedSet) const;

  /**
   *  @brief  Phase 6 - write results back to Pandora:
   *          AddToCluster for seeded hits, Cluster::Create for unseeded.
   */
  pandora::StatusCode WriteBack(std::vector<SeededCluster>& seeds, std::vector<UnseededCluster>& unseeded) const;

  // ================================================================== //
  //  Helper functions
  // ================================================================== //

  /**
   *  @brief  Helper: project (theta, phi) direction onto the HCAL inner face.
   *          Uses barrel vs. endcap geometry based on polar angle.
   */
  void ProjectThetaPhiToHcalFace(float theta, float phi, pandora::CartesianVector& out) const;

  /// Opening-angle distance: 1 - cos(angle between cluster direction and hit position)
  float OpeningAngleDist(const pandora::CartesianVector& clusDir, const pandora::CartesianVector& hitPos) const;

  /// Recompute the DR-corrected energy and log-weighted barycenter of any
  /// ClusterState (or derived type: SeededCluster, UnseededCluster) from its hits.
  void ComputeClusterState(ClusterState& cs) const;

  bool IsHcalHit(const pandora::CaloHit* h) const;
  bool IsCherenkovHit(const pandora::CaloHit* h) const;

  /**
   *  @brief  Retrieve the HCAL inner surface (barrel inner radius, endcap inner |z|) from the
   *          Pandora geometry.  Called from ReadSettings: the client creates the geometry before
   *          PandoraApi::ReadSettings, so the sub-detectors are already registered by then.
   */
  pandora::StatusCode ReadHcalGeometry();

  // ================================================================== //
  //  Configurable parameters
  // ================================================================== //
  std::string m_outputClusterListName = "EcalSeededClusters"; ///< Output cluster list name
  bool m_useDualReadout = true;                               ///< - use dual-readout info for seed selection and growth
  /// Take the helix magnetic field at the interaction point instead of at the seed's track state at
  /// the calorimeter.  Default true: the current field description is exactly zero outside the
  /// solenoid volume.  Switch off once a realistic (non-uniform) field map is available.
  bool m_useBFieldAtIP = true;
  float m_pThres = 1.0f;         ///< GeV - min track p for charged seed
  float m_hardThreshold = 0.00f; ///< GeV - min HCAL hit energy for pool
  float m_growThreshold = 0.00f; ///< GeV - min hit energy that propagates BFS

  // HCAL inner surface: retrieved from the Pandora geometry in ReadSettings (not configurable).
  float m_hcalSurfaceR = 0.f;     ///< mm - HCAL barrel inner radius (HCAL_BARREL inner R)
  float m_hcalSurfaceZHalf = 0.f; ///< mm - HCAL endcap inner |z|    (HCAL_ENDCAP inner Z)
  bool m_hasHcalEndcap = true;    ///< false for a barrel-only geometry (no HCAL_ENDCAP registered)

  float m_maxSearchDeltaR = 0.3f;  ///< rad - per-seed vicinity radius
  float m_seedCircleRadius = 50.f; ///< mm - initial seed circle radius for seeding BFS growth
  // Strip-anchored seed-circle search windows (FULL widths in rad; the cut is +-0.5x).
  // Neutral seeds project from the shower-biased ECAL centroid (offset to the true HCAL
  // energy is larger and phi-dominated), so they use a wider window than charged seeds,
  // which project from the precise track.
  float m_stripDThetaCharged = 0.10f; ///< rad - full theta window, charged seeds
  float m_stripDPhiCharged = 0.20f;   ///< rad - full phi   window, charged seeds
  float m_stripDThetaNeutral = 0.20f; ///< rad - full theta window, neutral seeds
  float m_stripDPhiNeutral = 0.20f;   ///< rad - full phi   window, neutral seeds
  float m_adjacencyRadius = 25.f;     ///< mm - BFS adjacency radius
  float m_rhoSeed = 0.08f;            ///< GeV - min scint sum in seed circle
  float m_rhoUnseeded = 0.08f;        ///< GeV - min scint sum for unseeded candidate
  float m_w0 = 4.6f;                  ///< - log-weight cutoff for barycenter
  float m_chiHcal = 0.31f;            ///< - dual-readout chi for HCAL
  float m_chargedEcalOverPMax = 0.7f; ///< - skip HCAL growth for a charged seed whose E_DR_ECAL/p exceeds this
}; // class EcalSeededClusteringAlgorithm

} // namespace lc_content

#endif // ECAL_SEEDED_CLUSTERING_ALGORITHM_H
