/**
 *  @file   LCContent/src/MLInference/SatelliteAssignmentOnnxAlgorithm.cc
 *
 *  @brief  ONNX cluster-grouping satellite assignment (SPLIT policy).  Token set, feature layout
 *          and pfo-indexing reproduce the model's training tokenisation.
 */
#include "Api/PandoraContentApi.h"
#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/Track.h"
#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/VectorHelper.h"
#include "MLInference/SatelliteAssignmentOnnxAlgorithm.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

using namespace pandora;

namespace lc_content {

//------------------------------------------------------------------------------------------------------------------------------------------

SatelliteAssignmentOnnxAlgorithm::~SatelliteAssignmentOnnxAlgorithm() = default;

void SatelliteAssignmentOnnxAlgorithm::LoadModel() {
  m_session.reset(new OnnxSession(m_modelPath)); // one intra-op thread; invalid on a bad path
}

//------------------------------------------------------------------------------------------------------------------------------------------

SatelliteAssignmentOnnxAlgorithm::Hit SatelliteAssignmentOnnxAlgorithm::ExtractHit(const CaloHit* const pHit) {
  Hit hit;
  float r = 0.f, phi = 0.f, theta = 0.f;
  pHit->GetPositionVector().GetSphericalCoordinates(r, phi, theta); // a calo hit is never at the origin
  hit.theta = theta;
  hit.phi = phi;
  hit.energy = pHit->GetInputEnergy();
  hit.isEcal = (pHit->GetHadronicEnergy() <= 0.f); // ECAL hits carry no hadronic energy
  hit.isCher = (pHit->GetHitType() == DRC_CHEREN);
  hit.depth = hit.isEcal ? static_cast<float>(pHit->GetLayer()) : SENT; // ECAL layer index; SENT for HCAL
  return hit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SatelliteAssignmentOnnxAlgorithm::BuildInfo(const Cluster* const pCluster, ClusterInfo& info) const {
  CaloHitList caloHits;
  pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
  if (caloHits.empty())
    return false;

  double sumX = 0., sumY = 0., sumZ = 0., sumE = 0.;
  for (const CaloHit* const pHit : caloHits) {
    const float energy = pHit->GetInputEnergy();
    if (energy <= 0.f)
      continue;
    info.hits.push_back(ExtractHit(pHit));
    const CartesianVector& position = pHit->GetPositionVector();
    sumX += position.GetX() * energy;
    sumY += position.GetY() * energy;
    sumZ += position.GetZ() * energy;
    sumE += energy;
  }
  if (info.hits.empty() || sumE <= 0.)
    return false;

  const CartesianVector centroid(static_cast<float>(sumX / sumE), static_cast<float>(sumY / sumE),
                                 static_cast<float>(sumZ / sumE));
  if (centroid.GetMagnitude() <= std::numeric_limits<float>::epsilon())
    return false;

  info.pCluster = pCluster;
  info.centroidDir = centroid.GetUnitVector();

  // any track reaching the calorimeter makes this a charged primary (track state taken there, no
  // projection; no |p| cut -- matches the training tokenisation); the highest-momentum track is used
  const Track* pBestTrack = nullptr;
  float bestMomentum = 0.f;
  for (const Track* const pTrack : pCluster->GetAssociatedTrackList()) {
    if (!pTrack->ReachesCalorimeter())
      continue;
    const float momentum = pTrack->GetTrackStateAtCalorimeter().GetMomentum().GetMagnitude();
    if (momentum > bestMomentum) {
      bestMomentum = momentum;
      pBestTrack = pTrack;
    }
  }
  if (pBestTrack) {
    float r = 0.f, phi = 0.f, theta = 0.f;
    pBestTrack->GetTrackStateAtCalorimeter().GetPosition().GetSphericalCoordinates(r, phi, theta);
    info.hasTrack = true;
    info.trkTheta = theta;
    info.trkPhi = phi;
    info.trkP = bestMomentum;
  }
  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// The ONE place the model runs.  Emits the candidate's own hits (pfo 0) plus, for each in-cone
// neighbour cluster, its top-K hits (and a track token if it is tracked) under a shared pfo index.
// The ONNX graph returns one affinity logit per neighbour pfo; we sigmoid them and return one
// NeighbourAff per neighbour (its cluster, whether it is tracked, and the affinity).

bool SatelliteAssignmentOnnxAlgorithm::ComputeNeighbourAffinities(const std::vector<ClusterInfo>& clusters,
                                                                  const int candIndex,
                                                                  std::vector<NeighbourAff>& out) const {
  out.clear();
  if (!m_session || !m_session->IsValid())
    return false;

  // Reference frame = candidate's energy-weighted centroid direction; all angles are relative to it.
  const ClusterInfo& candidate = clusters[candIndex];
  float r = 0.f, anchorPhi = 0.f, anchorTheta = 0.f;
  candidate.centroidDir.GetSphericalCoordinates(r, anchorPhi, anchorTheta);

  std::vector<float> tokenData;
  std::vector<int64_t> typeIds, pfoIds;

  // A calo-hit token: [dTheta, dPhi, depth, ln E, isCher, isEcal].
  auto addHit = [&](const Hit& hit, const int type, const int pfo) {
    tokenData.push_back(hit.theta - anchorTheta);
    tokenData.push_back(VectorHelper::deltaPhi(hit.phi, anchorPhi));
    tokenData.push_back(hit.depth);
    tokenData.push_back(std::log(std::max(hit.energy, 1e-9f)));
    tokenData.push_back(hit.isCher ? 1.f : 0.f);
    tokenData.push_back(hit.isEcal ? 1.f : 0.f);
    typeIds.push_back(type);
    pfoIds.push_back(pfo);
  };

  // A track token: the ln|p| carrier at the calorimeter; depth/cher/ecal are sentinels.
  auto addTrack = [&](const ClusterInfo& cluster, const int pfo) {
    tokenData.push_back(cluster.trkTheta - anchorTheta);
    tokenData.push_back(VectorHelper::deltaPhi(cluster.trkPhi, anchorPhi));
    tokenData.push_back(SENT);
    tokenData.push_back(std::log(std::max(cluster.trkP, 1e-6f)));
    tokenData.push_back(SENT);
    tokenData.push_back(SENT);
    typeIds.push_back(TYPE_TRACK);
    pfoIds.push_back(pfo);
  };

  // Keep only the top-K hits by energy (dynamic sequence length -> no padding).
  auto topByEnergy = [](std::vector<Hit> hits, const int k) {
    if (static_cast<int>(hits.size()) > k) {
      std::nth_element(hits.begin(), hits.begin() + k, hits.end(),
                       [](const Hit& a, const Hit& b) { return a.energy > b.energy; });
      hits.resize(static_cast<std::size_t>(k));
    }
    return hits;
  };

  // Own tokens live under pfo 0.
  for (const Hit& hit : topByEnergy(candidate.hits, m_budgetOwn))
    addHit(hit, TYPE_OWN, 0);

  // Every in-cone neighbour (no cap on their number) gets a 1-based pfo index; the k-th neighbour
  // maps to output slot k-1.
  std::vector<std::pair<const Cluster*, bool>> neighbourOfPfo; // per pfo: (cluster, hasTrack)
  int pfo = 0;
  for (int j = 0; j < static_cast<int>(clusters.size()); ++j) {
    if (j == candIndex)
      continue;
    if (candidate.centroidDir.GetCosOpeningAngle(clusters[j].centroidDir) <= m_coneCos)
      continue;

    ++pfo;
    const ClusterInfo& neighbour = clusters[j];
    for (const Hit& hit : topByEnergy(neighbour.hits, m_budgetNbr))
      addHit(hit, TYPE_OTHER, pfo);
    if (neighbour.hasTrack)
      addTrack(neighbour, pfo);

    neighbourOfPfo.emplace_back(neighbour.pCluster, neighbour.hasTrack);
  }

  const int nNeighbours = pfo;
  if (nNeighbours == 0) // nothing in cone -> no affinities
    return false;

  // Pack the inputs (order matches the model's inputs) and run via the shared ONNX session.
  const int64_t nTokens = static_cast<int64_t>(typeIds.size());
  const std::vector<OnnxSession::Input> inputs = {
      {{1, nTokens, N_FEAT}, OnnxSession::Input::Type::Float, tokenData.data()},
      {{1, nTokens}, OnnxSession::Input::Type::Int64, typeIds.data()},
      {{1, nTokens}, OnnxSession::Input::Type::Int64, pfoIds.data()},
  };
  const std::vector<std::vector<float>> outputs = m_session->Run(inputs);
  if (outputs.empty() || outputs[0].size() < static_cast<std::size_t>(nNeighbours))
    return false; // inference failed / short output

  const std::vector<float>& logits = outputs[0]; // length == nNeighbours (dynamic)
  out.reserve(nNeighbours);
  for (int k = 0; k < nNeighbours; ++k) // pfo k+1 -> output slot k
    out.push_back({neighbourOfPfo[k].first, neighbourOfPfo[k].second, 1.f / (1.f + std::exp(-logits[k]))}); // sigmoid
  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Snapshot every cluster in the current list into a ClusterInfo (hits, centroid, track state).
// Clusters BuildInfo rejects (no usable hits / degenerate centroid) are dropped.

pandora::StatusCode SatelliteAssignmentOnnxAlgorithm::BuildClusterCache(std::vector<ClusterInfo>& clusters) const {
  const ClusterList* pClusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

  for (const Cluster* const pCluster : *pClusterList) {
    ClusterInfo info;
    if (this->BuildInfo(pCluster, info))
      clusters.push_back(std::move(info));
  }
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// The merge candidates are the clusters whose charged/neutral fate this algorithm decides: available,
// trackless, and not already flagged a photon.  A cluster with a calo-reaching track is a charged
// primary (a target, never a candidate); photons are handled by the upstream photon-ID chain.

void SatelliteAssignmentOnnxAlgorithm::CollectCandidates(const std::vector<ClusterInfo>& clusters,
                                                         std::vector<int>& candidates) const {
  for (int i = 0; i < static_cast<int>(clusters.size()); ++i) {
    const Cluster* const pCluster = clusters[i].pCluster;

    const bool isCandidate = pCluster->IsAvailable() && !clusters[i].hasTrack && (PHOTON != pCluster->GetParticleId());
    if (isCandidate)
      candidates.push_back(i);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// SPLIT policy (unchanged, frozen): each candidate is merged into the TRACKED neighbour it most wants
// to join, provided that affinity clears AffinityThreshold.  Trackless neighbours are ignored here
// (they are the neutral pass's business).  Nothing is applied yet -- we only record child -> parent.

void SatelliteAssignmentOnnxAlgorithm::PlanChargedMerges(
    const std::vector<ClusterInfo>& clusters, const std::vector<int>& candidates,
    const std::map<int, std::vector<NeighbourAff>>& affinities,
    std::map<const Cluster*, const Cluster*>& chargedMerges) const {
  for (const int i : candidates) {
    const auto iter = affinities.find(i);
    if (iter == affinities.end()) // inference failed / no neighbour -> keep neutral
      continue;

    // Highest-affinity tracked neighbour.
    const Cluster* pBestParent = nullptr;
    float bestAffinity = -1.f;
    for (const NeighbourAff& n : iter->second) {
      if (n.hasTrack && n.affinity > bestAffinity) {
        bestAffinity = n.affinity;
        pBestParent = n.pCluster;
      }
    }

    const Cluster* const pCandidate = clusters[i].pCluster;
    if (bestAffinity > m_affinityThreshold && pBestParent && pBestParent != pCandidate)
      chargedMerges[pCandidate] = pBestParent;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Union-find over neutral<->neutral affinities (transitive: A~B, B~C => one object).  Only clusters
// that stayed neutral participate (those that went charged are frozen and excluded), and only neutral
// neighbours are ever unioned -- a charged cluster is never a target.  Each resulting component is
// emitted with its biggest cluster first, so ApplyNeutralGroups merges the rest into it.

void SatelliteAssignmentOnnxAlgorithm::PlanNeutralGroups(
    const std::vector<ClusterInfo>& clusters, const std::vector<int>& candidates,
    const std::map<int, std::vector<NeighbourAff>>& affinities, const std::set<const Cluster*>& becameCharged,
    std::vector<std::vector<const Cluster*>>& neutralGroups) const {
  // Seed the disjoint-set with the surviving neutral clusters (each its own parent).
  std::map<const Cluster*, const Cluster*> parent;
  for (const int i : candidates) {
    const Cluster* const pCluster = clusters[i].pCluster;
    if (!becameCharged.count(pCluster))
      parent[pCluster] = pCluster;
  }

  // Disjoint-set find with path-halving (only ever called on keys present in `parent`).
  auto find = [&parent](const Cluster* x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };

  // Union each surviving candidate with its neutral neighbours whose affinity clears the threshold.
  for (const int i : candidates) {
    const Cluster* const pCandidate = clusters[i].pCluster;
    if (!parent.count(pCandidate)) // candidate went charged -> not in the set
      continue;

    const auto iter = affinities.find(i);
    if (iter == affinities.end())
      continue;

    for (const NeighbourAff& n : iter->second) {
      const bool isNeutralSurvivor = !n.hasTrack && parent.count(n.pCluster);
      if (isNeutralSurvivor && n.affinity > m_neutralAffinityThreshold)
        parent[find(pCandidate)] = find(n.pCluster);
    }
  }

  // Collect components; keep only real groups (>= 2 members) and put the biggest cluster first.
  std::map<const Cluster*, std::vector<const Cluster*>> components;
  for (const auto& entry : parent)
    components[find(entry.first)].push_back(entry.first);

  const auto energy = [](const Cluster* c) { return c->GetHadronicEnergy() + c->GetElectromagneticEnergy(); };
  for (auto& entry : components) {
    std::vector<const Cluster*>& members = entry.second;
    if (members.size() < 2)
      continue;

    std::size_t repIdx = 0;
    for (std::size_t k = 1; k < members.size(); ++k)
      if (energy(members[k]) > energy(members[repIdx]))
        repIdx = k;

    std::swap(members[0], members[repIdx]); // representative = biggest, merged INTO
    neutralGroups.push_back(members);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Execute the planned charged merges.  The parent is a tracked cluster (never itself a merge child),
// so it cannot be stale here; a null/unavailable parent is skipped defensively.

pandora::StatusCode
SatelliteAssignmentOnnxAlgorithm::ApplyChargedMerges(const std::map<const Cluster*, const Cluster*>& merges) const {
  for (const auto& childParent : merges) {
    const Cluster* const pParent = childParent.second;
    if (!pParent || !pParent->IsAvailable())
      continue;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::MergeAndDeleteClusters(*this, pParent, childParent.first));
  }
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Execute the planned neutral groups: merge every non-representative member into the representative
// (members[0]).  Groups are disjoint, so no member is ever touched twice and the representative is
// never deleted.

pandora::StatusCode
SatelliteAssignmentOnnxAlgorithm::ApplyNeutralGroups(const std::vector<std::vector<const Cluster*>>& groups) const {
  for (const std::vector<const Cluster*>& members : groups) {
    const Cluster* const pRep = members.front();

    for (const Cluster* const pChild : members) {
      if (pChild == pRep)
        continue;

      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pRep, pChild));
    }
  }
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode SatelliteAssignmentOnnxAlgorithm::Run() {
  if (!m_session || !m_session->IsValid())
    return STATUS_CODE_FAILURE;

  // Snapshot the clusters and pick the trackless non-photon candidates.
  std::vector<ClusterInfo> clusters;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildClusterCache(clusters));

  std::vector<int> candidates;
  this->CollectCandidates(clusters, candidates);

  // Score each candidate ONCE: its affinity to every in-cone neighbour (charged and neutral alike).
  // Both passes below read these same scores -- no second inference (matches the offline validation).
  std::map<int, std::vector<NeighbourAff>> affinities;
  for (const int i : candidates) {
    std::vector<NeighbourAff> neighbours;
    if (this->ComputeNeighbourAffinities(clusters, i, neighbours))
      affinities.emplace(i, std::move(neighbours));
  }

  // Pass 1: plan the (frozen) charged merges; remember which candidates thereby became charged.
  std::map<const Cluster*, const Cluster*> chargedMerges;
  this->PlanChargedMerges(clusters, candidates, affinities, chargedMerges);

  std::set<const Cluster*> becameCharged;
  for (const auto& childParent : chargedMerges)
    becameCharged.insert(childParent.first);

  // Pass 2: plan the neutral<->neutral grouping among the clusters that stayed neutral.
  std::vector<std::vector<const Cluster*>> neutralGroups;
  if (m_groupNeutral)
    this->PlanNeutralGroups(clusters, candidates, affinities, becameCharged, neutralGroups);

  // Apply charged first (freezes the charged assignment), then neutral among the survivors.
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ApplyChargedMerges(chargedMerges));
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ApplyNeutralGroups(neutralGroups));

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode SatelliteAssignmentOnnxAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                           XmlHelper::ReadValue(xmlHandle, "ModelPath", m_modelPath)); // required

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "AffinityThreshold", m_affinityThreshold));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ConeHalfAngle", m_coneHalfAngle));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetOwn", m_budgetOwn));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetNeighbour", m_budgetNbr));

  // 2nd pass: neutral<->neutral union-find grouping (defaults on / 0.72; see Grouping_Satellite_Findings.md)
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "GroupNeutral", m_groupNeutral));
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "NeutralAffinityThreshold", m_neutralAffinityThreshold));

  m_coneCos = std::cos(m_coneHalfAngle);
  this->LoadModel();
  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
