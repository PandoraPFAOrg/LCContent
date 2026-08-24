#include "LCPfoConstruction/IdeaPfoCreationAlgorithm.h"
#include "LCHelpers/SortingHelper.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include "LCHelpers/VectorHelper.h"

#include <algorithm>     // std::min / std::max
#include <cmath>         // std::cos
#include <limits>        // std::numeric_limits
#include <unordered_map> // regional cache

namespace lc_content {

IdeaPfoCreationAlgorithm::IdeaPfoCreationAlgorithm() {}

const pandora::Track* IdeaPfoCreationAlgorithm::FindBestAssociatedTrack(const pandora::Cluster* aClus) const {
  const auto& associatedTrackList = aClus->GetAssociatedTrackList();

  if (associatedTrackList.empty())
    return nullptr;

  auto sortByTrackClusterDistance = [&aClus](const pandora::Track* a, const pandora::Track* b) {
    return SortingHelper::SortTracksByDistance(a, b, aClus);
  };

  std::vector<const pandora::Track*> tracksVector;
  tracksVector.insert(tracksVector.end(), associatedTrackList.begin(), associatedTrackList.end());

  std::sort(tracksVector.begin(), tracksVector.end(), sortByTrackClusterDistance);

  return tracksVector.front();
} // FindBestAssociatedTrack

float IdeaPfoCreationAlgorithm::EstimateScintEnergy(const pandora::Cluster* aClus) const {
  // E_S = ecalS + hcalS: scintillation sum with the same accumulation as the dual-readout plugin
  // (ECAL scint at EM scale, HCAL scint at hadronic scale), isolated hits included.
  pandora::CaloHitList hits;
  aClus->GetOrderedCaloHitList().FillCaloHitList(hits);
  const auto& isolated = aClus->GetIsolatedCaloHitList();

  float eS = 0.f;
  for (const pandora::CaloHitList& hitList : {hits, isolated}) {
    for (const auto* hit : hitList) {
      if (hit->GetHitType() != pandora::DRC_SCINT)
        continue;
      const float em = hit->GetElectromagneticEnergy(); // > 0 for ECAL
      eS += (em > 0.f) ? em : hit->GetHadronicEnergy(); // else HCAL scint (hadronic scale)
    }
  }

  return eS;
} // EstimateScintEnergy

pandora::StatusCode IdeaPfoCreationAlgorithm::GetDualReadoutEnergy(const pandora::Cluster* aClus, float& energy) const {
  const auto* pEnergyCorrections = PandoraContentApi::GetPlugins(*this)->GetEnergyCorrections();

  if (!pEnergyCorrections) {
    std::cerr << "IdeaPfoCreationAlgorithm: Unable to retrieve EnergyCorrections plugin!" << std::endl;
    return pandora::STATUS_CODE_FAILURE;
  }

  float correctedEm = 0.f, correctedHad = 0.f;
  pEnergyCorrections->MakeEnergyCorrections(aClus, correctedEm, correctedHad);
  energy = correctedHad; // use hadronic (dual-readout corrected) energy

  return pandora::STATUS_CODE_SUCCESS;
} // GetDualReadoutEnergy

pandora::StatusCode IdeaPfoCreationAlgorithm::IsEmShower(const pandora::Cluster* aClus, bool& isEmShower) const {
  const auto* pParticleId = PandoraContentApi::GetPlugins(*this)->GetParticleId();

  if (!pParticleId) {
    std::cerr << "IdeaPfoCreationAlgorithm: Unable to retrieve ParticleId plugin!" << std::endl;
    return pandora::STATUS_CODE_FAILURE;
  }

  isEmShower = pParticleId->IsEmShower(aClus);

  return pandora::STATUS_CODE_SUCCESS;
} // IsEmShower

pandora::StatusCode IdeaPfoCreationAlgorithm::Run() {
  const pandora::PfoList* pPfoList = nullptr;
  std::string pfoListName;

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                           PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

  // for now rely on the "Calo-driven" way
  pandora::ClusterList chargedClusters, neutralClusters;

  const pandora::ClusterList* clusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, clusterList));

  for (auto iter = clusterList->begin(); iter != clusterList->end(); ++iter) {
    const pandora::Cluster* const aClus = *iter;
    const auto& associatedTrackList = aClus->GetAssociatedTrackList();

    if (!associatedTrackList.empty()) {
      chargedClusters.push_back(aClus);
    } else {
      neutralClusters.push_back(aClus);
    }
  }

  // retrieve PID plugins
  // const auto* pParticleId = PandoraContentApi::GetPlugins(*this)->GetParticleId();

  // if (!pParticleId) {
  //   std::cerr << "IdeaPfoCreationAlgorithm: Unable to retrieve ParticleId plugin!" << std::endl;
  //   return pandora::STATUS_CODE_FAILURE;
  // }

  // fetch the track list up front (needed by CreatePfoFromTrack)
  const pandora::TrackList* trackList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, trackList));

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateElectronCandidates(chargedClusters));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateChargedHadronCandidates(chargedClusters));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreatePhotonCandidates(neutralClusters));
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateNeutralHadronCandidates(neutralClusters));

  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreatePfoFromTrack(trackList));

  if (!pPfoList->empty()) {
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::SaveList<pandora::ParticleFlowObject>(*this, m_outputPfoListName));
    PANDORA_RETURN_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::ReplaceCurrentList<pandora::ParticleFlowObject>(*this, m_outputPfoListName));
  }

  return pandora::STATUS_CODE_SUCCESS;
} // Run

float IdeaPfoCreationAlgorithm::CaloSigma(float p, bool isEm) const {
  if (p <= 0.f)
    return std::numeric_limits<float>::max();

  const float stoch = isEm ? m_stochasticEm : m_stochasticHad;
  const float cst = isEm ? m_constantEm : m_constantHad;
  return std::max(m_sigmaFloor, stoch / std::sqrt(p) + cst);
}

pandora::StatusCode
IdeaPfoCreationAlgorithm::CreateElectronCandidates(const pandora::ClusterList& /*clusterList*/) const {
  // No electron ID yet -> all track-associated clusters are handled by the charged-hadron branch
  // (electrons inside jets are rare).  Disabled (not removed): original body commented out below.
  return pandora::STATUS_CODE_SUCCESS;
  /*
  // loop over clusters
  for (auto iter = clusterList.begin(); iter != clusterList.end(); ++ iter) {
    const pandora::Cluster* const aClus = *iter;

    bool isEmShower = aClus->GetParticleId() == pandora::PHOTON; // use photon flag set upstream for now
    // PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->IsEmShower(aClus, isEmShower));

    if (!isEmShower)
      continue; // TODO electron ID can be different from the photon ID

    const auto& associatedTrackList = aClus->GetAssociatedTrackList();
    bool hasMutlipleTracks = associatedTrackList.size() > 1;

    // create PFO
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_clusterList.push_back(aClus);

    for (const auto* aTrack : associatedTrackList) {
      // add all associated tracks for now
      pfoParameters.m_trackList.push_back(aTrack);
    }

    const auto* bestTrack = FindBestAssociatedTrack(aClus);

    // TODO apply proper electron ID and use different energy estimation for electrons
    // use cluster property if there are multiple tracks
    float energy = bestTrack->GetMomentumAtDca().GetMagnitude();
    float clusterEnergy = 0.f;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->GetDualReadoutEnergy(aClus, clusterEnergy));
    auto momentum = bestTrack->GetMomentumAtDca();
    const float eop = clusterEnergy / energy;
    const float sigmaEm = CaloSigma(energy, isEmShower);
    bool failEoverP = eop > 1.f + m_nSigma * sigmaEm;

    if (hasMutlipleTracks || failEoverP) {
      // fall back to cluster energy if there are multiple tracks or fail E/p cut
      energy = clusterEnergy;

      // calculate energy and momentum
      auto clusterPos = aClus->GetCentroid(aClus->GetInnerPseudoLayer()); // TODO use beamspot
      auto unitVec = clusterPos.GetUnitVector();
      momentum = unitVec * energy;
    }

    pfoParameters.m_energy = energy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::E_MINUS);
    pfoParameters.m_charge = bestTrack->GetCharge();
    pfoParameters.m_particleId = (pfoParameters.m_charge.Get() > 0) ?
  pandora::PdgTable::GetParticlePdgCode(pandora::E_PLUS) : pandora::PdgTable::GetParticlePdgCode(pandora::E_MINUS);

    // TODO add vertex

    // Create the pfo
    const pandora::ParticleFlowObject* aPFO = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this,
  pfoParameters, aPFO)); } // cluster loop

  return pandora::STATUS_CODE_SUCCESS;
  */
} // CreateElectronCandidates

pandora::StatusCode
IdeaPfoCreationAlgorithm::CreateChargedHadronCandidates(const pandora::ClusterList& clusterList) const {
  // ---- regional-excess significance gate for the charged-hadron energy -------------------------------
  //   Single-track cluster: over a 3D opening-angle cone (m_coneOpeningAngle), sum E_DR of ALL clusters
  //   (charged + neutral) and p of all associated tracks.  Form the regional neutral excess and its
  //   significance against the calo resolution:
  //       S = (sumEdr - sumP) / ( sigma_c(sumEdr) * sumEdr )
  //   S > m_significanceK : a real neutral is merged into the region -> keep the precise track p on the
  //                         charged PFO and EMIT the FULL cluster excess (E_DR - p) as a separate
  //                         clusterless neutral-hadron PFO along the cluster centroid (a genuine particle).
  //   otherwise           : pure charged -> charged PFO = track p, no neutral.
  //   The charged energy is NEVER inflated (the tracker measurement is left untouched); the recovered
  //   neutral is always a distinct object with its physical calo-minus-track energy.
  //   Multi-track cluster -> whole E_DR at the centroid (the per-track gate is undefined with >1 track).
  //   Regional cache (E_DR, centroid, sum track p per cluster) is built ONCE so the cone loop never
  //   re-loops the calo hits.  The track (IP) direction for the charged PFO is what wins the downstream
  //   jet assignment (curling low-pT tracks: the momentum tangent at the calo points to the wrong jet).

  const pandora::ClusterList* fullClusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, fullClusterList));

  struct RegInfo {
    pandora::CartesianVector dir;
    float edr;
    float sumP;
  };
  std::unordered_map<const pandora::Cluster*, RegInfo> cache;
  cache.reserve(fullClusterList->size());
  for (const auto* c : *fullClusterList) {
    float edr = 0.f;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->GetDualReadoutEnergy(c, edr));
    float sumP = 0.f;
    for (const auto* trk : c->GetAssociatedTrackList())
      sumP += trk->GetMomentumAtDca().GetMagnitude(); // sum ALL associated tracks (not best)
    cache.emplace(c, RegInfo{c->GetCentroid(c->GetInnerPseudoLayer()).GetUnitVector(), edr, sumP});
  }

  const float cosCone = std::cos(m_coneOpeningAngle);

  // loop over the (track-associated) charged clusters and build the PFOs
  for (auto iter = clusterList.begin(); iter != clusterList.end(); ++iter) {
    const pandora::Cluster* const aClus = *iter;

    const auto& associatedTrackList = aClus->GetAssociatedTrackList();
    const bool hasMutlipleTracks = associatedTrackList.size() > 1;

    // create PFO
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_clusterList.push_back(aClus);

    for (const auto* aTrack : associatedTrackList) {
      // add all associated tracks for now
      pfoParameters.m_trackList.push_back(aTrack);
    }

    const auto* bestTrack = FindBestAssociatedTrack(aClus);

    const RegInfo& self = cache.at(aClus);
    const float clusterEnergy = self.edr; // this cluster's E_DR (from the cache)
    const float trackP = bestTrack->GetMomentumAtDca().GetMagnitude();

    // regional E/p over the opening-angle cone (self included: cos = 1 > cosCone)
    float sumEdr = 0.f, sumP = 0.f;
    for (const auto& kv : cache) {
      if (self.dir.GetCosOpeningAngle(kv.second.dir) > cosCone) {
        sumEdr += kv.second.edr;
        sumP += kv.second.sumP;
      }
    }
    const float sigmaReg = CaloSigma(sumEdr, /*isEm=*/false); // fractional hadronic sigma at the regional E_DR
    const float regExcess = sumEdr - sumP;                    // regional neutral energy (calo above tracks)
    const bool gateFires = regExcess > m_significanceK * sigmaReg * sumEdr; // S = regExcess/(sigma_c*sumEdr) > k
    const float clusterExcess = clusterEnergy - trackP;                     // this cluster's calo-minus-track energy

    float energy = trackP;                         // charged PFO keeps the exact track p (never inflated)
    auto momentum = bestTrack->GetMomentumAtDca(); // track (IP) momentum, magnitude p
    float newNeutralEn = 0.f;

    if (hasMutlipleTracks) {
      energy = clusterEnergy; // multi-track -> whole E_DR at centroid
      momentum = self.dir * energy;
    } else if (gateFires && clusterExcess > 0.f) {
      newNeutralEn = clusterExcess; // gate fires: emit the FULL excess as a real neutral
    }
    // else: pure charged -> energy = trackP, no neutral

    pfoParameters.m_energy = energy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass =
        pandora::PdgTable::GetParticleMass(pandora::PI_PLUS); // TODO pion mass for now (revisit after K-pi separation)
    pfoParameters.m_charge = bestTrack->GetCharge();
    pfoParameters.m_particleId = (pfoParameters.m_charge.Get() > 0)
                                     ? pandora::PdgTable::GetParticlePdgCode(pandora::PI_PLUS)
                                     : pandora::PdgTable::GetParticlePdgCode(pandora::PI_MINUS);

    // TODO add vertex

    // Create the charged pfo
    const pandora::ParticleFlowObject* aPFO = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, aPFO));

    // Emitted co-axial neutral (gate fired) -> separate neutral-hadron PFO along the cluster centroid.
    // It carries NO cluster: Pandora enforces exclusive cluster->PFO ownership (the charged PFO already
    // owns this cluster -- sharing it throws STATUS_CODE_NOT_ALLOWED), and this neutral is an INFERRED
    // energy from the charged cluster's calo excess with no calo of its own.  Jet reco uses its
    // 4-momentum (energy = excess, centroid direction).
    if (newNeutralEn > 0.f) {
      PandoraContentApi::ParticleFlowObject::Parameters nhParams;
      nhParams.m_energy = newNeutralEn;
      nhParams.m_momentum = self.dir * newNeutralEn;
      nhParams.m_mass = pandora::PdgTable::GetParticleMass(pandora::K_LONG);
      nhParams.m_charge = 0.f;
      nhParams.m_particleId = pandora::PdgTable::GetParticlePdgCode(pandora::K_LONG);
      const pandora::ParticleFlowObject* nhPFO = nullptr;
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::ParticleFlowObject::Create(*this, nhParams, nhPFO));
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
} // CreateChargedHadronCandidates

pandora::StatusCode IdeaPfoCreationAlgorithm::CreatePhotonCandidates(const pandora::ClusterList& clusterList) const {
  // loop over clusters
  for (auto iter = clusterList.begin(); iter != clusterList.end(); ++iter) {
    const pandora::Cluster* const aClus = *iter;

    // use photon flag set upstream
    if (aClus->GetParticleId() != pandora::PHOTON)
      continue;

    // Photon energy = E_S (scintillation only).  For EM showers S ~= C, so the dual-readout
    // combination E_DR = (S - chi*C)/(1 - chi) only adds Cherenkov sampling noise; the plain
    // scintillation sum E_S = ecalS + hcalS is the lower-variance photon estimator.
    const float photonEnergy = this->EstimateScintEnergy(aClus);

    auto clusterPos = aClus->GetCentroid(aClus->GetInnerPseudoLayer()); // TODO use beamspot
    auto unitVec = clusterPos.GetUnitVector();
    auto photonMomentum = unitVec * photonEnergy;

    // create PFO
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_clusterList.push_back(aClus);
    pfoParameters.m_energy = photonEnergy;
    pfoParameters.m_momentum = photonMomentum;
    pfoParameters.m_mass = 0.f;
    pfoParameters.m_charge = 0.f;
    pfoParameters.m_particleId = pandora::PdgTable::GetParticlePdgCode(pandora::PHOTON);

    // Create the pfo
    const pandora::ParticleFlowObject* aPFO = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, aPFO));
  }

  return pandora::STATUS_CODE_SUCCESS;
} // CreatePhotonCandidates

pandora::StatusCode
IdeaPfoCreationAlgorithm::CreateNeutralHadronCandidates(const pandora::ClusterList& clusterList) const {
  // loop over clusters
  for (auto iter = clusterList.begin(); iter != clusterList.end(); ++iter) {
    const pandora::Cluster* const aClus = *iter;

    if (aClus->GetParticleId() == pandora::PHOTON)
      continue; // flagged as photon: skip

    // run dual-readout correction for energy estimation
    float clusterEnergy = 0.f;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->GetDualReadoutEnergy(aClus, clusterEnergy));

    // Flat hadronic-response calibration on the E_DR neutral-hadron energy.  E_DR under-responds for
    // these (mostly untracked-charged-as-neutral) hadron clusters by ~10% at the event level; the
    // jet-energy-budget / NH-scale-scan study shows a constant 1.10 lift de-biases the di-jet mass to
    // ~91 GeV with NO mass-resolution cost.  (The recovered-neutral K_LONG PFOs in the charged branch
    // are regression-calibrated to the true scale, so they are NOT scaled.)  XML "NeutralHadEnergyScale".
    clusterEnergy *= m_neutralHadScale;

    // create PFO
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_clusterList.push_back(aClus);

    // calculate energy and momentum
    auto clusterPos = aClus->GetCentroid(aClus->GetInnerPseudoLayer()); // TODO use beamspot
    auto unitVec = clusterPos.GetUnitVector();
    auto momentum = unitVec * clusterEnergy;

    pfoParameters.m_energy = clusterEnergy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::K_LONG);
    pfoParameters.m_charge = 0.f;
    pfoParameters.m_particleId = aClus->GetParticleId(); // pandora::PdgTable::GetParticlePdgCode(pandora::K_LONG);

    // Create the pfo
    const pandora::ParticleFlowObject* aPFO = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, aPFO));
  }

  return pandora::STATUS_CODE_SUCCESS;
} // CreateNeutralHadronCandidates

pandora::StatusCode IdeaPfoCreationAlgorithm::CreatePfoFromTrack(const pandora::TrackList* /*trackList*/) const {
  // DISABLED: clusterless track-only ("looper") PFO emission removed -- it DOUBLE-COUNTS energy.
  // These low-pT tracks have no associated cluster, but truth shows ~91% of them DO reach the calo and
  // deposit into CLUSTERED cells (already counted as NH/photon); the track is then added again here as a
  // muon PFO (~0.8 GeV/evt double-count).  Root cause is a track-cluster ASSOCIATION failure: Pandora's
  // single-helix extrapolation gives a (0,0,0) AtCalorimeter state for 58% of these curlers, and lands
  // >15 mm off the nearest ECAL hit for most of the rest, so the ECAL 10 mm match cannot fire.  The
  // (rejected) forward-|eta| hypothesis was 0%.  Proper fix is upstream: a real track finder + Kalman fit
  // with material/extrapolation, then these clusters associate and become charged.  Until then, do NOT
  // emit them.  (The co-axial-neutral K_LONG from the E_DR/p gate in the charged branch is unaffected.)
  return pandora::STATUS_CODE_SUCCESS;
  /*
  const pandora::ClusterList* clusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::GetCurrentList(*this, clusterList));

  // loop over tracks and create a PFO for each track not already used in the cluster-based PFOs
  for (const auto* track : *trackList) {
    if (track->HasAssociatedCluster())
      continue;

    float pt, phi, pz;
    track->GetMomentumAtDca().GetCylindricalCoordinates(pt, phi, pz);

    if (pt > m_ptCut)
      continue; // it should have reached the calorimeter

    // create PFO
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_trackList.push_back(track);

    // calculate energy and momentum
    float energy = track->GetMomentumAtDca().GetMagnitude();
    auto momentum = track->GetMomentumAtDca();

    pfoParameters.m_energy = energy;
    pfoParameters.m_momentum = momentum;
    pfoParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::MU_MINUS); // muon for now
    pfoParameters.m_charge = track->GetCharge();
    pfoParameters.m_particleId = (pfoParameters.m_charge.Get() > 0) ?
  pandora::PdgTable::GetParticlePdgCode(pandora::MU_PLUS) : pandora::PdgTable::GetParticlePdgCode(pandora::MU_MINUS);

    // Create the pfo
    const pandora::ParticleFlowObject* aPFO = nullptr;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this,
  pfoParameters, aPFO));
  }

  return pandora::STATUS_CODE_SUCCESS;
  */
} // CreatePfoFromTrack

pandora::StatusCode IdeaPfoCreationAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "NSigma", m_nSigma));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "StochasticTermHad", m_stochasticHad));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "ConstantTermHad", m_constantHad));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "StochasticTermEm", m_stochasticEm));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "ConstantTermEm", m_constantEm));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "PtCut", m_ptCut));

  // regional-excess significance gate on the charged-hadron energy
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "ConeOpeningAngle", m_coneOpeningAngle));

  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "SignificanceK", m_significanceK));

  // flat neutral-hadron energy calibration (independent of the regression; default 1.10)
  PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
                                  pandora::XmlHelper::ReadValue(xmlHandle, "NeutralHadEnergyScale", m_neutralHadScale));

  return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lc_content
