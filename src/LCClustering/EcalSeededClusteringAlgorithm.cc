/**
 *  @file   LCContent/src/LCClustering/EcalSeededClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the ECAL-seeded HCAL clustering algorithm.
 *
 *  Phase summary
 *  -------------
 *   1a  SelectSeeds    - filter ECAL clusters into SeededCluster list
 *   1b  BuildHcalPool    - cache all HCAL hits above HardThreshold
 *   2a  ExtrapolateToHcal  - project each seed onto the HCAL inner face
 *   2b  BuildVicinities  - restrict per-seed work to dR < MaxSearchDeltaR
 *   3   FormSeedCircles +  - initial seed-circle frontier + contest resolution
 *       ResolveSeedCircleContests
 *   4   GrowSeeds     - layer-by-layer BFS within each seed's vicinity
 *   5   FormUnseededClusters - BFS+union-find for hits not claimed in Phase 4
 *   6   WriteBack      - AddToCluster / Cluster::Create
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Api/PandoraContentApi.h"
#include "Objects/Track.h"
#include "Objects/Cluster.h"
#include "Objects/CaloHit.h"
#include "Objects/Helix.h"

#include "LCClustering/EcalSeededClusteringAlgorithm.h"
#include "LCHelpers/VectorHelper.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

namespace lc_content {

// ====================================================================== //
//  Run
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::Run() {
  // ---- get Pandora input lists ----
  // assume all existing clusters are ECAL clusters
  const pandora::ClusterList *pEcalClusterList  = nullptr;
  const pandora::CaloHitList *pCurrentCaloHitList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::GetCurrentList(*this, pEcalClusterList));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::GetCurrentList(*this, pCurrentCaloHitList));

  // ---- Phase 1a: ECAL seed selection ----
  std::vector<SeededCluster> seeds;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      SelectSeeds(*pEcalClusterList, seeds));

  // ---- Phase 1b: HCAL hit pool ----
  std::vector<HitCache> hcalPool; // this is the owner of the cached hit info
  BuildHcalPool(*pCurrentCaloHitList, hcalPool);

  if (hcalPool.empty()) // nothing to do
    return pandora::STATUS_CODE_SUCCESS;

  // ---- Phases 2–4: seed-driven work (skipped when no seeds) ----
  pandora::CaloHitSet assignedSet;

  if (!seeds.empty()) {
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
        FormSeededClusters(seeds, hcalPool, assignedSet));
  }

  // ---- Phase 5: unseeded clusters from leftover hits ----
  // Runs regardless of whether any ECAL seeds were present.
  std::vector<UnseededCluster> unseeded;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      FormUnseededClusters(unseeded, hcalPool, assignedSet));

  // ---- Phase 6: write back to Pandora ----
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      WriteBack(seeds, unseeded));

  return pandora::STATUS_CODE_SUCCESS;
}

// ====================================================================== //
//  Phase 1a - SelectSeeds
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::SelectSeeds(const pandora::ClusterList &ecalClusters,
                                                               std::vector<SeededCluster> &seeds) const {
  // loop over ECAL clusters and filter out
  for (const pandora::Cluster *const pClus : ecalClusters) {
    if (!pClus->IsAvailable())
      continue;

    SeededCluster s;
    s.pEcalCluster = pClus;

    const auto &tracks = pClus->GetAssociatedTrackList();
    if (tracks.empty()) {
      // Only grow HCAL behind neutral clusters tagged as neutral hadrons by the
      // external PID (cluster particleId == K_LONG, set upstream from the edm4hep
      // cluster type). Photons, satellites and unidentified neutral clusters are
      // skipped.
      if (pClus->GetParticleId() != pandora::K_LONG)
        continue;

      s.isCharged  = false;
      s.pBestTrack = nullptr;
    } else {
      // Charged cluster: take the highest-momentum associated track
      const pandora::Track *bestTrack = nullptr;
      float bestP = 0.f;
      for (const auto *t : tracks) {
        const float p = t->GetMomentumAtDca().GetMagnitude();

        if (p > bestP) {
          bestP = p;
          bestTrack = t;
        }
      } // loop over associated tracks

      // Drop seeds below momentum threshold
      if (bestP < m_pThres)
        continue;

      // E_DR_ECAL/p gate: a charged cluster that already deposited most of its
      // momentum in the ECAL has little real HCAL behind it, so growing would only
      // attach a neighbour's HCAL (fake). The DR (dual-readout) correction is
      // returned in the hadronic slot; this algorithm runs before HCAL growing so
      // the cluster is ECAL-only and corrHad is just E_DR_ECAL.
      float eEcal = 0.f;
      if (m_useDualReadout) {
        const pandora::EnergyCorrections *const pCorr =
            PandoraContentApi::GetPlugins(*this)->GetEnergyCorrections();
        if (pCorr == nullptr)
          return pandora::STATUS_CODE_FAILURE;

        float corrEm = 0.f, corrHad = 0.f;
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
            pCorr->MakeEnergyCorrections(pClus, corrEm, corrHad));
        eEcal = corrHad;
      } else {
        eEcal = pClus->GetElectromagneticEnergy();
      }

      if (eEcal >= m_chargedEcalOverPMax * bestP)
        continue;

      s.isCharged  = true;
      s.pBestTrack = bestTrack;
    }

    seeds.push_back(std::move(s));
  } // loop over ECAL clusters

  return pandora::STATUS_CODE_SUCCESS;
} // SelectSeeds

// ====================================================================== //
//  Phase 1b - BuildHcalPool
// ====================================================================== //

void EcalSeededClusteringAlgorithm::BuildHcalPool(const pandora::CaloHitList &allHits,
                                                  std::vector<HitCache> &hcalPool) const {
  // loop over all hits and cache HCAL hits above threshold
  hcalPool.reserve(allHits.size());
  for (const pandora::CaloHit *const h : allHits) {
    if (h->GetHadronicEnergy() <= 0.f)
      continue; // skip ECAL
    if (h->GetHadronicEnergy() < m_hardThreshold)
      continue;

    HitCache hc;
    hc.pHit  = h;
    hc.pos   = h->GetPositionVector();
    hc.energy  = h->GetHadronicEnergy();

    float r;
    hc.pos.GetSphericalCoordinates(r, hc.phi, hc.theta);
    hc.isCheren = IsCherenkovHit(h);
    hcalPool.push_back(std::move(hc));
  }
} // BuildHcalPool

// ====================================================================== //
//  Phase 2a - ExtrapolateToHcalSurface
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::ExtrapolateToHcalSurface(SeededCluster &s) const {
  pandora::CartesianVector impact(0.f, 0.f, 0.f);

  if (s.isCharged && s.pBestTrack != nullptr) {
    // ---------------------------------------------------------------- //
    //  Helix extrapolation from track state at ECAL surface using
    //  pandora::Helix::GetPointOnCircle  (barrel)  and
    //  pandora::Helix::GetPointInZ       (endcap).
    // ---------------------------------------------------------------- //
    const pandora::TrackState& tsAtCalo = s.pBestTrack->GetTrackStateAtCalorimeter(); 
    const pandora::CartesianVector &rp  = tsAtCalo.GetPosition();
    const pandora::CartesianVector &mom = tsAtCalo.GetMomentum();
    const float pz = mom.GetZ();

    // Magnetic field for the helix.  UseBfieldAtIP (default) takes the field at the IP
    // since exactly zero magnetic field outside the solenoid can invalidate pandora::Helix.
    // With a realistic field map, switch the option off to evaluate the field
    // at the extrapolation START point (the track state at the calorimeter) instead.
    // Note that pandora::Helix itself assumes a constant field
    // a fair approximation over the short ECAL -> HCAL surface step.
    const pandora::CartesianVector bFieldPosition =
        m_useBFieldAtIP ? pandora::CartesianVector(0.f, 0.f, 0.f) : rp;
    const float bField = PandoraContentApi::GetPlugins(*this)->GetBFieldPlugin()->GetBField(bFieldPosition);

    // extrapolation based on the track state at calo
    const pandora::Helix helix(rp, mom, static_cast<float>(s.pBestTrack->GetCharge()), bField);

    // Try barrel intersection
    pandora::CartesianVector barrelPoint(0.f, 0.f, 0.f);
    const pandora::StatusCode scBarrel =
        helix.GetPointOnCircle(m_hcalSurfaceR, rp, barrelPoint);
    const float barrelTime = (barrelPoint.GetZ() - rp.GetZ()) / pz;

    // Try endcap intersection (correct sign from momentum z-component).  Skipped for a
    // barrel-only geometry: m_hcalSurfaceZHalf is the infinite sentinel there, and asking the
    // helix to intersect a plane at that z is not meaningful.
    bool endcapOk = false;
    pandora::CartesianVector endcapPoint(0.f, 0.f, 0.f);
    float endcapTime = 0.f;

    if (m_hasHcalEndcap) {
      const float targetZ = (pz >= 0.f) ? m_hcalSurfaceZHalf : -m_hcalSurfaceZHalf;
      const pandora::StatusCode scEndcap = helix.GetPointInZ(targetZ, rp, endcapPoint);
      endcapTime = (endcapPoint.GetZ() - rp.GetZ()) / pz;
      const float endcapR2 = endcapPoint.GetX()*endcapPoint.GetX() +
                             endcapPoint.GetY()*endcapPoint.GetY();
      endcapOk = (scEndcap == pandora::STATUS_CODE_SUCCESS) &&
                 (endcapR2 <= m_hcalSurfaceR * m_hcalSurfaceR) &&
                 (endcapTime > 0.f); // require forward extrapolation in time
    }

    const bool barrelOk = (scBarrel == pandora::STATUS_CODE_SUCCESS) &&
                          (std::fabs(barrelPoint.GetZ()) <= m_hcalSurfaceZHalf) &&
                          (barrelTime > 0.f); // require forward extrapolation in time

    if (barrelOk && endcapOk) {
      // Pick whichever is closer to the reference point in time
      impact = (barrelTime <= endcapTime) ? barrelPoint : endcapPoint;
    } else if (barrelOk) {
      impact = barrelPoint;
    } else if (endcapOk) {
      impact = endcapPoint;
    } else {
      std::cerr << "EcalSeededClusteringAlgorithm::ExtrapolateToHcalSurface - "
                << "Failed to extrapolate charged seed to HCAL surface" << std::endl;
      return pandora::STATUS_CODE_FAILURE;
    } // if barrel vs. endcap intersection

    // reject if the extrapolated point is
    // unreasonably far from the track state at calo
    float thImpact, phImpact, thCalo, phCalo, _;
    impact.GetSphericalCoordinates(_, phImpact, thImpact);
    rp.GetSphericalCoordinates(_, phCalo, thCalo);
    if (VectorHelper::AngularDR(thImpact, phImpact, thCalo, phCalo) > m_maxSearchDeltaR)
      return pandora::STATUS_CODE_FAILURE; // fail silently (likely low pT tracks)    
  } else {
    // Neutral seed: project ECAL centroid radially through origin
    const pandora::CartesianVector &centroid =
        s.pEcalCluster->GetCentroid(s.pEcalCluster->GetInnerPseudoLayer());
    float th, ph, _;
    centroid.GetSphericalCoordinates(_, ph, th);
    ProjectThetaPhiToHcalFace(th, ph, impact);
  } // if charged vs. neutral seed

  float _;
  s.hcalImpactPoint = impact;
  impact.GetSphericalCoordinates(_, s.phiSeed, s.thetaSeed);

  return pandora::STATUS_CODE_SUCCESS;
} // ExtrapolateToHcalSurface

void EcalSeededClusteringAlgorithm::ProjectThetaPhiToHcalFace(float theta,
                                                              float phi,
                                                              pandora::CartesianVector &out) const {
  const float thetaTrans = std::atan2(m_hcalSurfaceR, m_hcalSurfaceZHalf);
  const bool  barrel   = (theta > thetaTrans) && (theta < (static_cast<float>(M_PI) - thetaTrans));

  if (barrel) {
    const float rho = m_hcalSurfaceR;
    // for a barrel point: x=r*cos(phi), y=r*sin(phi), z=r/tan(theta)
    const float z   = rho / std::tan(theta);
    out = pandora::CartesianVector(rho * std::cos(phi), rho * std::sin(phi), z);
  } else {
    const float z   = (std::cos(theta) >= 0.f) ? m_hcalSurfaceZHalf : -m_hcalSurfaceZHalf;
    // z and tan(theta) always share the same sign here (both >0 for +z endcap,
    // both <0 for -z endcap), so z*tan(theta) is always positive.
    const float rho = z * std::tan(theta);
    out = pandora::CartesianVector(rho * std::cos(phi), rho * std::sin(phi), z);
  }
} // ProjectThetaPhiToHcalFace

// ====================================================================== //
//  Phase 2b - BuildVicinities
// ====================================================================== //

void EcalSeededClusteringAlgorithm::BuildVicinities(std::vector<SeededCluster> &seeds,
                                                    const std::vector<HitCache> &hcalPool) const {
  const float maxDR2 = m_maxSearchDeltaR * m_maxSearchDeltaR;
  for (auto &s : seeds) {
    s.vicinity.clear();

    if (!s.alive)
      continue;

    for (const auto &hc : hcalPool) {
      if (VectorHelper::AngularDR2(s.thetaSeed, s.phiSeed, hc.theta, hc.phi) < maxDR2)
        s.vicinity.push_back(&hc);
    } // loop over hcalPool
  } // loop over seeds
} // BuildVicinities

// ====================================================================== //
//  Phase 3a - FormSeedCircles
// ====================================================================== //

void EcalSeededClusteringAlgorithm::FormSeedCircles(std::vector<SeededCluster> &seeds,
                                                    SeedCircleMap &circles) const {
  for (std::size_t i = 0; i < seeds.size(); ++i) {
    SeededCluster &s = seeds[i];
    if (!s.alive)
      continue;

    // Strip search window (half-widths) around the projected impact direction.
    const float halfDTheta = 0.5f * (s.isCharged ? m_stripDThetaCharged : m_stripDThetaNeutral);
    const float halfDPhi   = 0.5f * (s.isCharged ? m_stripDPhiCharged   : m_stripDPhiNeutral);

    // Collect candidate anchor hits inside the (theta, phi) strip, ordered by
    // opening-angle distance to the projection so the FIRST one passing the
    // density cut is automatically the closest qualifying anchor.
    std::vector<const HitCache *> candidates;
    for (const HitCache *hc : s.vicinity) {
      if (std::fabs(hc->theta - s.thetaSeed) > halfDTheta)
        continue;
      if (std::fabs(VectorHelper::deltaPhi(hc->phi, s.phiSeed)) > halfDPhi)
        continue;
      candidates.push_back(hc);
    }
    std::sort(candidates.begin(), candidates.end(),
              [&s](const HitCache *a, const HitCache *b) {
                return VectorHelper::AngularDR2(s.thetaSeed, s.phiSeed, a->theta, a->phi) <
                       VectorHelper::AngularDR2(s.thetaSeed, s.phiSeed, b->theta, b->phi);
              });

    // For each candidate (closest first) build its m_seedCircleRadius disk; the
    // first whose scintillation sum passes m_rhoSeed becomes the seed circle.
    std::vector<const HitCache *> circle;
    bool found = false;
    for (const HitCache *anchor : candidates) {
      const float r_anchor = anchor->pos.GetMagnitude();
      const float cosCut = r_anchor / std::sqrt(r_anchor * r_anchor + m_seedCircleRadius * m_seedCircleRadius);

      circle.clear();
      float scintSum = 0.f;
      for (const HitCache *hc : s.vicinity) {
        if (hc->energy < m_growThreshold)
          continue;
        if (anchor->pos.GetCosOpeningAngle(hc->pos) > cosCut) {
          circle.push_back(hc);
          if (!hc->isCheren)
            scintSum += hc->energy;
        }
      } // loop over vicinity hits for this anchor's disk

      if (scintSum >= m_rhoSeed) {
        found = true;
        break; // closest qualifying anchor wins
      }
    } // loop over strip candidates

    // Kill seeds with no qualifying anchor in the strip
    if (!found) {
      s.alive = false;
      continue;
    }

    circles[static_cast<int>(i)] = std::move(circle);
  } // loop over seeds
} // FormSeedCircles

// ====================================================================== //
//  Phase 3b - ResolveSeedCircleContests
// ====================================================================== //

void EcalSeededClusteringAlgorithm::ResolveSeedCircleContests(std::vector<SeededCluster> &seeds,
                                                              SeedCircleMap &circles) const {
  // Build hit -> {seed indices} map
  std::map<const HitCache *, std::vector<int>> contestants;
  for (const auto &[idx, circle] : circles)
    for (const HitCache *hc : circle)
      contestants[hc].push_back(idx);

  // Resolve each contested hit with opening angle to the seed impact point
  for (const auto &[hc, owners] : contestants) {
    if (owners.size() < 2)
      continue;

    // Pick winner: lowest opening-angle distance
    int   winner = owners[0];
    float bestDist = std::numeric_limits<float>::max();
    for (int idx : owners) {
      if (!seeds[idx].alive)
        continue;

      const float d = OpeningAngleDist(seeds[idx].hcalImpactPoint, hc->pos);
      if (d < bestDist) { bestDist = d; winner = idx; }
    } // loop over contestants for this hit

    // Remove from all losers
    for (int idx : owners) {
      if (idx == winner)
        continue;

      auto &cl = circles[idx];
      cl.erase(std::remove(cl.begin(), cl.end(), hc), cl.end());

      if (cl.empty())
        seeds[idx].alive = false;
    } // loop over contestants for this hit
  } // loop over contested hits
} // ResolveSeedCircleContests

// ====================================================================== //
//  ComputeClusterState
// ====================================================================== //

void EcalSeededClusteringAlgorithm::ComputeClusterState(ClusterState &gc) const {
  float scintSum = 0.f, cherenSum = 0.f, totE = 0.f;
  for (const HitCache *hc : gc.hits) {
    totE += hc->energy;
    if (hc->isCheren) cherenSum += hc->energy;
    else              scintSum  += hc->energy;
  }

  gc.energy = m_useDualReadout
      ? (scintSum - m_chiHcal * cherenSum) / (1.f - m_chiHcal)
      : totE;

  pandora::CartesianVector bary(0.f, 0.f, 0.f);
  float totW = 0.f;
  if (totE > 0.f) {
    for (const HitCache *hc : gc.hits) {
      const float w = std::max(0.f, m_w0 + std::log(hc->energy / totE));
      if (w > 0.f) {
        bary += hc->pos * w;
        totW += w;
      }
    } // loop over hits in cluster
  }

  gc.barycenter = (totW > 0.f) ? bary * (1.f / totW)
                                : pandora::CartesianVector(0.f, 0.f, 0.f);
} // ComputeClusterState

// ====================================================================== //
//  Shared BFS growth engine - Grow()
// ====================================================================== //

std::vector<EcalSeededClusteringAlgorithm::ClusterState>
EcalSeededClusteringAlgorithm::Grow(GrowStrategy strategy,
                                    const std::vector<ClusterState> &clusters,
                                    pandora::CaloHitSet &assignedSet) const {
  std::vector<ClusterState> result(clusters); // working copy
  const int n = static_cast<int>(result.size());
  if (n == 0)
    return {};

  const float adjCut = m_adjacencyRadius; // interpreted as physical distance for cone angle

  // Per-cluster frontier: hits added in the last layer.
  // Initialised from each cluster's seed-circle hits.
  std::map<int, std::vector<const HitCache *>> frontier;
  for (int i = 0; i < n; ++i)
    frontier[i] = result[i].hits;

  // Union-find bookkeeping (MergeClusters only).
  std::vector<int> rep(n);
  std::iota(rep.begin(), rep.end(), 0);
  auto findRep = [&rep](int x) -> int {
    while (rep[x] != x) { rep[x] = rep[rep[x]]; x = rep[x]; }
    return x;
  };
  auto mergeReps = [&rep, &findRep](int a, int b) {
    a = findRep(a); b = findRep(b);
    if (a != b) rep[b] = a;
  };

  // Shared lambda: accumulate candidates from a (frontier, vicinity) pair into
  // the candidates map, keyed by the given cluster index / representative id.
  auto collectCandidates = [&](std::map<const HitCache *, std::vector<int>> &candidates,
                               int id,
                               const std::vector<const HitCache *> &front,
                               const auto &vicinity) {
    for (const HitCache *fhc : front) {
      // Max opening angle for adjacency: atan2(adjCut, r_fhc)
      const float r_fhc  = fhc->pos.GetMagnitude();
      const float cosCut = r_fhc / std::sqrt(r_fhc * r_fhc + adjCut * adjCut);

      for (const HitCache *nb : vicinity) {
        if (assignedSet.count(nb->pHit))
          continue;
        if (nb->energy < m_growThreshold)
          continue;
        if (fhc->pos.GetCosOpeningAngle(nb->pos) > cosCut)
          candidates[nb].push_back(id);
      }
    }
  };

  // Shared lambda: deduplicate a sorted-or-unsorted claimers vector in-place.
  auto deduplicateClaimers = [](std::vector<int> &claimers) {
    std::sort(claimers.begin(), claimers.end());
    claimers.erase(std::unique(claimers.begin(), claimers.end()), claimers.end());
  };

  bool anyAdded = true;
  while (anyAdded) {
    anyAdded = false;

    // ---------------------------------------------------------------- //
    //  Collect candidates, deduplicate, assign, and build next frontier.
    //  ResolveByDist: each cluster searches its own vicinity independently.
    //  MergeClusters: build per-representative vicinity union, then search.
    // ---------------------------------------------------------------- //
    std::map<int, std::vector<const HitCache *>> nextFrontier;

    if (strategy == GrowStrategy::ResolveByDist) {
      std::map<const HitCache *, std::vector<int>> candidates;
      for (int i = 0; i < n; ++i)
        collectCandidates(candidates, i, frontier[i], result[i].vicinity);

      for (auto &[hc, claimers] : candidates) {
        deduplicateClaimers(claimers);

        int   winner = claimers[0];
        float bestDist = std::numeric_limits<float>::max();
        for (int idx : claimers) {
          const float d = OpeningAngleDist(result[idx].barycenter, hc->pos);
          if (d < bestDist) { bestDist = d; winner = idx; }
        }
        assignedSet.insert(hc->pHit);
        result[winner].hits.push_back(hc);
        nextFrontier[winner].push_back(hc);
        anyAdded = true;
      } // loop over candidates

      // recompute cluster states
      for (const auto &[i, _] : nextFrontier)
        ComputeClusterState(result[i]);
    } else { // MergeClusters: collapse frontier and vicinity into per-representative sets
      std::map<int, std::vector<const HitCache *>> repFront;
      std::map<int, std::set<const HitCache *>> repVicinity;
      for (int i = 0; i < n; ++i) {
        const int r = findRep(i);
        for (const HitCache *hc : frontier[i])
          repFront[r].push_back(hc);
        for (const HitCache *hc : result[i].vicinity)
          repVicinity[r].insert(hc);
      } // loop over clusters to build repFront and repVicinity

      std::map<const HitCache *, std::vector<int>> candidates;
      for (const auto &[r, fr] : repFront)
        collectCandidates(candidates, r, fr, repVicinity[r]);

      for (auto &[hc, claimers] : candidates) {
        deduplicateClaimers(claimers);

        int winner = findRep(claimers[0]);
        for (std::size_t j = 1; j < claimers.size(); ++j)
          mergeReps(winner, claimers[j]);
        const int ri = findRep(winner);
        assignedSet.insert(hc->pHit);
        result[ri].hits.push_back(hc);
        nextFrontier[ri].push_back(hc);
        anyAdded = true;
      } // loop over candidates
    } // if ResolveByDist vs. MergeClusters

    frontier = std::move(nextFrontier);
  } // while anyAdded

  // MergeClusters: fold non-representative clusters into their root
  // (hits + vicinity), then erase the now-empty entries.
  if (strategy == GrowStrategy::MergeClusters) {
    for (int i = 0; i < n; ++i) {
      const int r = findRep(i);
      if (r == i)
        continue;

      // Merge hits and vicinity into representative cluster
      for (const HitCache *hc : result[i].hits)
        result[r].hits.push_back(hc);
      for (const HitCache *hc : result[i].vicinity)
        result[r].vicinity.push_back(hc);

      // Clear the non-representative cluster's hits and vicinity.
      // Do NOT erase here — erasing inside the loop invalidates indices.
      result[i].hits.clear();
      result[i].vicinity.clear();
    } // loop over clusters to merge into representatives

    // Single erase pass after the loop: remove all now-empty non-representatives.
    result.erase(
        std::remove_if(result.begin(), result.end(),
                       [](const ClusterState &gc) { return gc.hits.empty(); }),
        result.end());
  } // if MergeClusters

  return result;
} // Grow

// ====================================================================== //
//  Phases 2–4 - FormSeededClusters
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::FormSeededClusters(std::vector<SeededCluster> &seeds,
                                                                      const std::vector<HitCache> &hcalPool,
                                                                      pandora::CaloHitSet &assignedSet) const {
  // ---- Phase 2a: project seeds onto HCAL surface ----
  for (auto &s : seeds) {
    if (ExtrapolateToHcalSurface(s) != pandora::STATUS_CODE_SUCCESS)
      s.alive = false;
  }

  // ---- Phase 2b: per-seed vicinity pools ----
  BuildVicinities(seeds, hcalPool);

  // ---- Phase 3: seed circles + contest ----
  SeedCircleMap circles;
  FormSeedCircles(seeds, circles);
  ResolveSeedCircleContests(seeds, circles);

  // ---- Phase 4: BFS growth via shared Grow() ----
  // Single pass: initialise ClusterStates from post-contest circles,
  // mark circle hits as assigned, and copy each seed's vicinity.
  std::vector<int> aliveIdx; // map initStates index -> seeds index
  std::vector<ClusterState> initStates;

  for (int i = 0; i < static_cast<int>(seeds.size()); ++i) {
    if (!seeds[i].alive)
      continue;
    aliveIdx.push_back(i);

    ClusterState cs;
    cs.hits     = circles[i];       // HitCache* directly - no map lookup needed
    cs.vicinity = seeds[i].vicinity; // search space for this seed
    for (const HitCache *hc : cs.hits)
      assignedSet.insert(hc->pHit);

    ComputeClusterState(cs);
    initStates.push_back(std::move(cs));
  } // loop over alive seeds

  const auto grown = Grow(GrowStrategy::ResolveByDist, initStates, assignedSet);

  // Write grown state back into each SeededCluster (inherits ClusterState).
  for (int j = 0; j < static_cast<int>(aliveIdx.size()); ++j) {
    auto &s = seeds[aliveIdx[j]];
    s.hits = grown[j].hits;
    ComputeClusterState(s);
  }

  return pandora::STATUS_CODE_SUCCESS;
} // FormSeededClusters

// ====================================================================== //
//  Phase 5 - FormUnseededClusters
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::FormUnseededClusters(std::vector<UnseededCluster> &unseeded,
                                                                        const std::vector<HitCache> &hcalPool,
                                                                        pandora::CaloHitSet &assignedSet) const {
  // Collect remaining HCAL hits not consumed by Phase 4.
  std::vector<const HitCache *> remaining;
  remaining.reserve(hcalPool.size());
  for (const auto &hc : hcalPool) {
    if (!assignedSet.count(hc.pHit))
      remaining.push_back(&hc);
  }

  if (remaining.empty()) // nothing to do
    return pandora::STATUS_CODE_SUCCESS;

  // Sort by descending energy for reproducible seeding order.
  std::sort(remaining.begin(), remaining.end(),
            [](const HitCache *a, const HitCache *b) { return a->energy > b->energy; });

  // ------------------------------------------------------------------ //
  //  Step A: find candidate seed centers - local scint circle sum >= rhoUnseeded.
  //  Circle hits are claimed greedily in energy order to avoid overlap.
  // ------------------------------------------------------------------ //
  pandora::CaloHitSet circleAssigned; // hits already claimed by a seed circle
  std::vector<ClusterState> initStates;
  const float maxDR2 = m_maxSearchDeltaR * m_maxSearchDeltaR;

  for (const HitCache *center : remaining) {
    // Max opening angles for this center hit: atan2(cut, r_center)
    const float r_center   = center->pos.GetMagnitude();
    const float cosSeedCut = r_center / std::sqrt(r_center * r_center + m_seedCircleRadius * m_seedCircleRadius);

    float scintSum = 0.f;
    ClusterState gc;
    for (const HitCache *nb : remaining) {
      const float dr2 = VectorHelper::AngularDR2(center->theta, center->phi, nb->theta, nb->phi);

      // Vicinity: all hits in range, regardless of threshold or circle ownership
      if (dr2 < maxDR2)
        gc.vicinity.push_back(nb);

      // Seed circle: above threshold and within seedCircleRadius, not already claimed
      if (circleAssigned.count(nb->pHit))
        continue;
      if (nb->energy < m_growThreshold)
        continue;

      const float cosA = center->pos.GetCosOpeningAngle(nb->pos);
      // form seed circle candidate
      if (cosA > cosSeedCut) {
        if (!nb->isCheren)
          scintSum += nb->energy;

        gc.hits.push_back(nb);
      }
    } // 2nd loop over remaining hits

    if (scintSum < m_rhoUnseeded)
      continue;

    for (const HitCache *hc : gc.hits)
      circleAssigned.insert(hc->pHit);

    ComputeClusterState(gc);
    initStates.push_back(std::move(gc));
  } // 1st loop over remaining hits

  if (initStates.empty())
    return pandora::STATUS_CODE_SUCCESS;

  // ------------------------------------------------------------------ //
  //  BFS with merge-on-contest via shared Grow().
  //  Pre-assign circle hits so Grow() won't re-claim them.
  // ------------------------------------------------------------------ //
  for (const auto &gc : initStates)
    for (const HitCache *hc : gc.hits)
      assignedSet.insert(hc->pHit);

  auto result = Grow(GrowStrategy::MergeClusters, initStates, assignedSet);

  // Emit one UnseededCluster per surviving root.
  for (auto &gc : result) {
    if (gc.hits.empty())
      continue;

    UnseededCluster u;
    u.hits = std::move(gc.hits);
    ComputeClusterState(u);
    unseeded.push_back(std::move(u));
  }

  return pandora::STATUS_CODE_SUCCESS;
} // FormUnseededClusters

// ====================================================================== //
//  Phase 6 - WriteBack
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::WriteBack(std::vector<SeededCluster> &seeds,
                                                             std::vector<UnseededCluster> &unseeded) const {
  // Extend each surviving ECAL cluster with its HCAL hits
  for (auto &s : seeds) {
    if (!s.alive || s.hits.empty())
      continue;

    for (const HitCache *hc : s.hits) {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
          PandoraContentApi::AddToCluster(*this, s.pEcalCluster, hc->pHit));
    }
  }

  // Get current list name before proceeding
  std::string inputClusterListName;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::GetCurrentListName<pandora::Cluster>(*this, inputClusterListName));

  // Recreate cluster list within the pandora framework
  // (this will allow creating new objects)
  const pandora::ClusterList* pClusterList = nullptr;
  std::string clusterListName = m_outputClusterListName + "Tmp";

  PANDORA_RETURN_RESULT_IF(
      pandora::STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

  // Create new clusters for unseeded HCAL clusters
  for (auto &u : unseeded) {
    pandora::CaloHitList hitList;
    for (const HitCache *hc : u.hits)
      hitList.push_back(hc->pHit);

    const pandora::Cluster *pNewCluster = nullptr;
    PandoraContentApi::Cluster::Parameters params;
    params.m_caloHitList = hitList;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::Cluster::Create(*this, params, pNewCluster));
  }

  // need these to store the pandora cluster list
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, // save the previous list to the new output (seeded clusters)
                            PandoraContentApi::SaveList<pandora::Cluster>(*this, inputClusterListName, m_outputClusterListName));

  if (!unseeded.empty()) { // only replace the current list if we actually created new clusters
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, // save the current list to the new output (unseeded clusters)
                              PandoraContentApi::SaveList<pandora::Cluster>(*this, m_outputClusterListName));
  }

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, // make the new output the current list
                            PandoraContentApi::ReplaceCurrentList<pandora::Cluster>(*this, m_outputClusterListName));

  return pandora::STATUS_CODE_SUCCESS;
}

// ====================================================================== //
//  Helper functions
// ====================================================================== //

float EcalSeededClusteringAlgorithm::OpeningAngleDist(const pandora::CartesianVector &clusDir,
                                                      const pandora::CartesianVector &hitPos) const {
  const float r_hit  = hitPos.GetMagnitude();
  const float r_clus = clusDir.GetMagnitude();
  if (r_hit <= 0.f || r_clus <= 0.f)
    return std::numeric_limits<float>::max();

  float cosA = clusDir.GetCosOpeningAngle(hitPos);
  cosA = std::max(-1.f, std::min(1.f, cosA));
  return 1.f - cosA;
} // OpeningAngleDist

bool EcalSeededClusteringAlgorithm::IsHcalHit(const pandora::CaloHit *h) const {
  return h->GetHadronicEnergy() > 0.f;
}

bool EcalSeededClusteringAlgorithm::IsCherenkovHit(const pandora::CaloHit *h) const {
  return h->GetHitType() == pandora::DRC_CHEREN;
}

// ====================================================================== //
//  ReadHcalGeometry
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::ReadHcalGeometry() {
  const pandora::GeometryManager *const pGeometry = PandoraContentApi::GetGeometry(*this);

  if (!pGeometry)
    return pandora::STATUS_CODE_NOT_INITIALIZED;

  // GetSubDetector throws StatusCodeException(NOT_FOUND) if the sub-detector was never registered;
  // the algorithm manager catches it around ReadSettings and reports it.
  m_hcalSurfaceR = pGeometry->GetSubDetector(pandora::HCAL_BARREL).GetInnerRCoordinate();

  // A barrel-only geometry registers no HCAL endcap.  That is a valid configuration, not an error:
  // the z half-length is then effectively infinite, so every direction lands on the barrel.
  try {
    m_hcalSurfaceZHalf = pGeometry->GetSubDetector(pandora::HCAL_ENDCAP).GetInnerZCoordinate();
  } catch (const pandora::StatusCodeException &) {
    m_hasHcalEndcap    = false;
    m_hcalSurfaceZHalf = std::numeric_limits<float>::max();
  }

  if ((m_hcalSurfaceR <= 0.f) || (m_hcalSurfaceZHalf <= 0.f))
    return pandora::STATUS_CODE_INVALID_PARAMETER;

  return pandora::STATUS_CODE_SUCCESS;
}

// ====================================================================== //
//  ReadSettings
// ====================================================================== //

pandora::StatusCode EcalSeededClusteringAlgorithm::ReadSettings(
    const pandora::TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "UseDualReadout", m_useDualReadout));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "MomentumThreshold", m_pThres));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "HardThreshold", m_hardThreshold));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "GrowThreshold", m_growThreshold));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "MaxSearchDR", m_maxSearchDeltaR));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "SeedCircleR", m_seedCircleRadius));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "StripDThetaCharged", m_stripDThetaCharged));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "StripDPhiCharged", m_stripDPhiCharged));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "StripDThetaNeutral", m_stripDThetaNeutral));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "StripDPhiNeutral", m_stripDPhiNeutral));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "AdjacencyR", m_adjacencyRadius));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "SeedThreshold", m_rhoSeed));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "UnseededThreshold", m_rhoUnseeded));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "W0", m_w0));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "ChiHcal", m_chiHcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "ChargedEcalOverPMax", m_chargedEcalOverPMax));
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "UseBfieldAtIP", m_useBFieldAtIP));

  // The HCAL inner surface comes from the detector geometry, not from XML.
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ReadHcalGeometry());

  if (m_useDualReadout) {
    if (m_chiHcal >= 1.f || m_chiHcal <= 0.f) {
      std::cerr << "EcalSeededClusteringAlgorithm::ReadSettings: ChiHcal must be 0 < chi < 1";
      std::cerr << " to avoid singularity in dual-readout correction." << std::endl;

      return pandora::STATUS_CODE_INVALID_PARAMETER;
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lc_content
