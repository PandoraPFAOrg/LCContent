/**
 *  @file   LCContent/src/MLInference/ClusterNeutralPidAlgorithm.cc
 *
 *  @brief  Implementation of the neutral-PID (sat/NH/photon) ONNX tagging algorithm.
 *
 *  Port of k4RecCalorimeter ClusterNeutralPidOnnx.  Tokenisation, feature layout and
 *  decision rule are kept 1:1 with that algorithm and with photonId_hitcollect_multi.py.
 *  Pandora accessors give exact feature parity (no cellID decode needed):
 *    ln E      <- CaloHit::GetInputEnergy()        (the digi-channel energy)
 *    is_cher   <- CaloHit::GetHitType()==DRC_CHEREN
 *    depth     <- CaloHit::GetLayer()              (= the cellID layer/"depth" field)
 *    ECAL/HCAL <- CaloHit::GetHadronicEnergy()<=0  => ECAL
 *    dTheta,dPhi <- CaloHit::GetPositionVector()
 *    track     <- Track::GetTrackStateAtCalorimeter().GetPosition()/.GetMomentum()
 *
 *  The model output is already a per-class score (no softmax applied here); the photon
 *  decision is a direct score cut at PhotonThreshold.
 */
#include "MLInference/ClusterNeutralPidAlgorithm.h"

#include "LCHelpers/VectorHelper.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace pandora;

namespace lc_content {

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterNeutralPidAlgorithm::ClusterNeutralPidAlgorithm()
    : m_modelPath(""), m_w0(4.6f), m_coneLateralMm(250.f), m_budgetOwnEcal(200), m_budgetOtherEcal(150),
      m_budgetHcal(150), m_budgetTrack(20), m_photonThreshold(-1.f), m_maxTokens(0) {}

ClusterNeutralPidAlgorithm::~ClusterNeutralPidAlgorithm() = default;

//------------------------------------------------------------------------------------------------------------------------------------------

int ClusterNeutralPidAlgorithm::Decide(float s0, float s1, float s2) const {
  // raw model scores in order [satellite, neutral-hadron, photon]; photon is a score cut
  if (m_photonThreshold >= 0.f) {
    if (s2 > m_photonThreshold)
      return PHOTON;
    return (s1 > s0) ? K_LONG : UNKNOWN_PARTICLE_TYPE;
  }
  // no threshold configured: fall back to argmax over the three scores
  if (s2 >= s1 && s2 >= s0)
    return PHOTON;
  return (s1 > s0) ? K_LONG : UNKNOWN_PARTICLE_TYPE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterNeutralPidAlgorithm::Run() {
  const ClusterList* pClusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
  if (pClusterList->empty())
    return STATUS_CODE_SUCCESS;

  // ------------------------------------------------------------------
  // 1. Hit pool: the current CaloHitList is the superset of all hits.
  // ------------------------------------------------------------------
  const CaloHitList* pCaloHitList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

  std::vector<Hit> hits;
  hits.reserve(pCaloHitList->size());
  for (const CaloHit* const pHit : *pCaloHitList) {
    const float e = pHit->GetInputEnergy();
    if (e <= 0.f)
      continue;
    const CartesianVector& p = pHit->GetPositionVector();
    float r, phi, theta;
    p.GetSphericalCoordinates(r, phi, theta); // a calo hit is never at the origin
    const CartesianVector u = p.GetUnitVector();
    Hit hh;
    hh.ux = u.GetX();
    hh.uy = u.GetY();
    hh.uz = u.GetZ();
    hh.th = theta;
    hh.ph = phi;
    hh.e = e;
    hh.isEcal = (pHit->GetHadronicEnergy() <= 0.f);
    hh.cherFeat = (pHit->GetHitType() == DRC_CHEREN) ? 1.f : 0.f;
    hh.depthFeat = static_cast<float>(pHit->GetLayer()); // used only for ECAL tokens
    hh.pCaloHit = pHit;
    hits.push_back(hh);
  }

  // ------------------------------------------------------------------
  // 2. Track impacts at the calorimeter.
  // ------------------------------------------------------------------
  const TrackList* pTrackList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pTrackList));

  std::vector<TrackImpact> impacts;
  for (const Track* const pTrack : *pTrackList) {
    // skip before GetSphericalCoordinates: the calo track-state position is (0,0,0)
    // when the track does not reach the calorimeter (which would otherwise throw).
    if (!pTrack->ReachesCalorimeter())
      continue;
    const TrackState& ts = pTrack->GetTrackStateAtCalorimeter();
    const CartesianVector& pos = ts.GetPosition();
    float r, phi, theta;
    pos.GetSphericalCoordinates(r, phi, theta);
    const CartesianVector u = pos.GetUnitVector();
    impacts.push_back({u.GetX(), u.GetY(), u.GetZ(), theta, phi, ts.GetMomentum().GetMagnitude()});
  }

  // ------------------------------------------------------------------
  // 3. Per-cluster tokenisation + inference + PID write-out.
  // ------------------------------------------------------------------
  std::vector<float> tokens(static_cast<std::size_t>(m_maxTokens) * N_FEAT);
  std::vector<int64_t> typeIds(static_cast<std::size_t>(m_maxTokens));
  std::unique_ptr<bool[]> mask(new bool[m_maxTokens]);

  struct Row {
    std::array<float, N_FEAT> f;
    float intensity;
    int64_t typeId;
  };
  std::vector<Row> ownList, otherList, hcalList, trkList;

  for (const Cluster* const pCluster : *pClusterList) {
    // skip clusters already associated with a track (charged): their identity is
    // resolved by track-cluster matching, and the model excludes charged primaries.
    if (!pCluster->GetAssociatedTrackList().empty())
      continue;

    int particleId = UNKNOWN_PARTICLE_TYPE; // satellite / untokenisable -> UNKNOWN

    // ---- own-hit set + total energy (for the log-weighted barycenter) ----
    std::set<const CaloHit*> ownSet;
    float eTotal = 0.f;
    for (const auto& layerEntry : pCluster->GetOrderedCaloHitList()) {
      for (const CaloHit* const pHit : *layerEntry.second) {
        ownSet.insert(pHit);
        const float e = pHit->GetInputEnergy();
        if (e > 0.f)
          eTotal += e;
      }
    }

    // ---- W0 log-weighted barycenter: w_i = max(0, m_w0 + ln(e_i / eTotal)) ----
    // reproduces the k4 cluster position the model trained against.
    double sx = 0., sy = 0., sz = 0., sw = 0.;
    if (eTotal > 0.f) {
      for (const auto& layerEntry : pCluster->GetOrderedCaloHitList()) {
        for (const CaloHit* const pHit : *layerEntry.second) {
          const float e = pHit->GetInputEnergy();
          if (e <= 0.f)
            continue;
          const float w = m_w0 + std::log(e / eTotal);
          if (w <= 0.f)
            continue;
          const CartesianVector& p = pHit->GetPositionVector();
          sx += p.GetX() * w;
          sy += p.GetY() * w;
          sz += p.GetZ() * w;
          sw += w;
        }
      }
    }

    if (sw > 0.) {
      const CartesianVector centroid(static_cast<float>(sx / sw), static_cast<float>(sy / sw),
                                     static_cast<float>(sz / sw));
      float rC, phC, thC;
      centroid.GetSphericalCoordinates(rC, phC, thC); // log-weighted barycenter, at calo radius
      const CartesianVector uc = centroid.GetUnitVector();
      const float ucx = uc.GetX(), ucy = uc.GetY(), ucz = uc.GetZ();
      const float cosAlpha = std::cos(std::atan(m_coneLateralMm / rC));

      ownList.clear();
      otherList.clear();
      hcalList.clear();
      trkList.clear();

      for (const Hit& h : hits) {
        const bool inOwn = ownSet.count(h.pCaloHit) > 0;
        const float cosGeom = ucx * h.ux + ucy * h.uy + ucz * h.uz;
        const float dth = h.th - thC;
        const float dph = VectorHelper::deltaPhi(h.ph, phC);
        const float lnE = std::log(h.e);

        if (h.isEcal && inOwn) {
          ownList.push_back({{dth, dph, h.depthFeat, lnE, h.cherFeat}, h.e, TID_OWN});
        } else if (h.isEcal && !inOwn && cosGeom >= cosAlpha) {
          otherList.push_back({{dth, dph, h.depthFeat, lnE, h.cherFeat}, h.e, TID_OTHER});
        } else if (!h.isEcal && cosGeom >= cosAlpha) {
          hcalList.push_back({{dth, dph, FEATURE_ABSENT, lnE, h.cherFeat}, h.e, TID_HCAL});
        }
      }

      for (const TrackImpact& imp : impacts) {
        const float cosGeom = ucx * imp.ux + ucy * imp.uy + ucz * imp.uz;
        if (cosGeom < cosAlpha)
          continue;
        const float dth = imp.th - thC;
        const float dph = VectorHelper::deltaPhi(imp.ph, phC);
        const float lnP = std::log(std::max(imp.p, 1e-6f));
        trkList.push_back({{dth, dph, FEATURE_ABSENT, lnP, FEATURE_ABSENT}, imp.p, TID_TRACK});
      }

      auto trim = [](std::vector<Row>& v, int k) {
        if (static_cast<int>(v.size()) <= k)
          return;
        std::nth_element(v.begin(), v.begin() + k, v.end(),
                         [](const Row& a, const Row& b) { return a.intensity > b.intensity; });
        v.resize(static_cast<std::size_t>(k));
      };
      trim(ownList, m_budgetOwnEcal);
      trim(otherList, m_budgetOtherEcal);
      trim(hcalList, m_budgetHcal);
      trim(trkList, m_budgetTrack);

      std::fill(tokens.begin(), tokens.end(), 0.f);
      std::fill(typeIds.begin(), typeIds.end(), int64_t{0});
      std::fill(mask.get(), mask.get() + m_maxTokens, true); // true = padding (masked)

      std::size_t idx = 0;
      auto emit = [&](const std::vector<Row>& v) {
        for (const Row& r : v) {
          const std::size_t off = idx * N_FEAT;
          for (int k = 0; k < N_FEAT; ++k)
            tokens[off + k] = r.f[k];
          typeIds[idx] = r.typeId;
          mask[idx] = false;
          ++idx;
        }
      };
      emit(ownList);
      emit(otherList);
      emit(hcalList);
      emit(trkList);

      const int64_t maxTok = static_cast<int64_t>(m_maxTokens);
      const std::vector<OnnxSession::Input> inputs = {
          {{1, maxTok, N_FEAT}, OnnxSession::Input::Type::Float, tokens.data()},
          {{1, maxTok}, OnnxSession::Input::Type::Int64, typeIds.data()},
          {{1, maxTok}, OnnxSession::Input::Type::Bool, mask.get()},
      };
      const std::vector<std::vector<float>> outputs = m_session->Run(inputs);
      if (outputs.empty() || outputs[0].size() < static_cast<std::size_t>(N_CLASSES)) {
        particleId = UNKNOWN_PARTICLE_TYPE; // inference failed -> UNKNOWN
      } else {
        const std::vector<float>& scores = outputs[0]; // [sat, NH, photon]
        particleId = Decide(scores[0], scores[1], scores[2]);
      }
    }

    // ---- write the verdict as the cluster ParticleId ----
    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = particleId;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterNeutralPidAlgorithm::LoadModel() {
  if (m_modelPath.empty()) {
    std::cout << "ClusterNeutralPidAlgorithm: ModelPath is empty." << std::endl;
    return STATUS_CODE_INVALID_PARAMETER;
  }
  m_session.reset(new OnnxSession(m_modelPath));
  if (!m_session->IsValid())
    return STATUS_CODE_FAILURE;

  if (m_session->InputCount() != 3 || m_session->OutputCount() != 1) {
    std::cout << "ClusterNeutralPidAlgorithm: expected 3 inputs and 1 output." << std::endl;
    return STATUS_CODE_FAILURE;
  }

  const std::vector<std::int64_t> tokensShape = m_session->InputShape(0);
  if (tokensShape.size() == 3 && tokensShape[1] > 0 && tokensShape[1] != static_cast<std::int64_t>(m_maxTokens)) {
    std::cout << "ClusterNeutralPidAlgorithm: model MAX_TOKENS=" << tokensShape[1] << " != configured budget sum "
              << m_maxTokens << "." << std::endl;
    return STATUS_CODE_FAILURE;
  }
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterNeutralPidAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelPath", m_modelPath));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "W0", m_w0));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ConeLateralMm", m_coneLateralMm));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetOwnEcal", m_budgetOwnEcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetOtherEcal", m_budgetOtherEcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetHcal", m_budgetHcal));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "BudgetTrack", m_budgetTrack));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "PhotonThreshold", m_photonThreshold));

  m_maxTokens = m_budgetOwnEcal + m_budgetOtherEcal + m_budgetHcal + m_budgetTrack;
  return this->LoadModel();
}

} // namespace lc_content
