/**
 *  @file   LCContent/src/LCParticleId/ForwardPhotonIdAlgorithm.cc
 *
 *  @brief  Implementation of the cut-based forward-photon tagging algorithm.
 */
#include "LCParticleId/ForwardPhotonIdAlgorithm.h"

#include "Pandora/AlgorithmHeaders.h"

#include <cmath>

using namespace pandora;

namespace lc_content {

//------------------------------------------------------------------------------------------------------------------------------------------

bool ForwardPhotonIdAlgorithm::CollectHcalHits(const Cluster *const pCluster,
                                               std::vector<HcalHit> &hits, float &eTotal) const {
  hits.clear();
  eTotal = 0.f;
  for (const auto &layerEntry : pCluster->GetOrderedCaloHitList()) {
    for (const CaloHit *const pHit : *layerEntry.second) {
      const float eHad = pHit->GetHadronicEnergy();
      if (eHad <= 0.f)                       // ECAL hit -> cluster is not HCAL-only
        return false;
      const CartesianVector &p = pHit->GetPositionVector();
      hits.push_back({p.GetX(), p.GetY(), p.GetZ(), eHad,
                      pHit->GetHitType() == DRC_CHEREN});
      eTotal += eHad;
    }
  }
  return (!hits.empty() && eTotal > 0.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ForwardPhotonIdAlgorithm::W0Centroid(const std::vector<HcalHit> &hits, const float eTotal,
                                          CartesianVector &centroid) const {
  double sx = 0., sy = 0., sz = 0., sw = 0.;
  for (const HcalHit &h : hits) {
    const float w = m_w0 + std::log(h.e / eTotal);
    if (w <= 0.f) continue;
    sx += h.x * w; sy += h.y * w; sz += h.z * w; sw += w;
  }
  if (sw <= 0.)
    return false;
  centroid = CartesianVector(static_cast<float>(sx / sw),
                             static_cast<float>(sy / sw),
                             static_cast<float>(sz / sw));
  return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ForwardPhotonIdAlgorithm::CherToScintRatio(const std::vector<HcalHit> &hits) const {
  float S = 0.f, C = 0.f;
  for (const HcalHit &h : hits) { if (h.isCher) C += h.e; else S += h.e; }
  if (S <= 0.f)
    return -1.f;                             // S == 0 -> skip (per design)
  return C / S;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ForwardPhotonIdAlgorithm::DRcorrHcal(const std::vector<HcalHit> &hits) const {
  float S = 0.f, C = 0.f;
  for (const HcalHit &h : hits) { if (h.isCher) C += h.e; else S += h.e; }
  return (S > 0.f && C > 0.f) ? (S - m_chiHcal * C) / (1.f - m_chiHcal) : S;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ForwardPhotonIdAlgorithm::CoreFraction(const std::vector<HcalHit> &hits,
                                             const CartesianVector &centroid,
                                             const float rCentroid) const {
  const float eDrTotal = this->DRcorrHcal(hits);
  if (eDrTotal <= 0.f)
    return 0.f;
  const CartesianVector u = centroid.GetUnitVector();
  const float ux = u.GetX(), uy = u.GetY(), uz = u.GetZ();
  const float cosCone = std::cos(std::atan(m_coreMm / rCentroid));
  std::vector<HcalHit> coreHits;
  for (const HcalHit &h : hits) {
    const float rh = std::sqrt(h.x*h.x + h.y*h.y + h.z*h.z);
    if (rh > 0.f && (ux*h.x + uy*h.y + uz*h.z) / rh > cosCone)
      coreHits.push_back(h);
  }
  return this->DRcorrHcal(coreHits) / eDrTotal;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ForwardPhotonIdAlgorithm::Run() {
  const ClusterList *pClusterList = nullptr;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
      PandoraContentApi::GetCurrentList(*this, pClusterList));

  for (const Cluster *const pCluster : *pClusterList) {
    if (!pCluster->IsAvailable())
      continue;
    if (!pCluster->GetAssociatedTrackList().empty())
      continue;                              // charged -> not a candidate

    const int particleId = pCluster->GetParticleId();
    if (PHOTON == particleId || K_LONG == particleId)
      continue;                              // already identified -> leave alone

    // ---- HCAL-only hit collection (also yields the total energy) ----
    std::vector<HcalHit> hcalHits;
    float eTotal = 0.f;
    if (!this->CollectHcalHits(pCluster, hcalHits, eTotal))
      continue;

    // ---- W0 log-weighted centroid ----
    CartesianVector centroid(0.f, 0.f, 0.f);
    if (!this->W0Centroid(hcalHits, eTotal, centroid))
      continue;
    float rC, phiC, thetaC;
    centroid.GetSphericalCoordinates(rC, phiC, thetaC);
    if (rC <= 0.f)
      continue;

    // ---- cut 0 (cheap, rejects most clusters): forward  |eta| > MinAbsEta ----
    const float th  = std::min(std::max(thetaC, 1e-6f),
                               static_cast<float>(M_PI) - 1e-6f);
    const float eta = -std::log(std::tan(0.5f * th));
    if (std::fabs(eta) <= m_minAbsEta)
      continue;

    // ---- cut 1: Cherenkov-richness  C/S > MinCS  (S == 0 -> skip) ----
    const float cs = this->CherToScintRatio(hcalHits);
    if (cs <= m_minCS)
      continue;

    // ---- cut 2: lateral compactness  coreFraction > MinCoreFraction ----
    if (this->CoreFraction(hcalHits, centroid, rC) <= m_minCoreFrac)
      continue;

    // ---- all cuts passed: tag as photon ----
    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = PHOTON;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ForwardPhotonIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "MinCS", m_minCS));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "MinAbsEta", m_minAbsEta));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "CoreMm", m_coreMm));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "MinCoreFraction", m_minCoreFrac));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "W0", m_w0));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ChiHcal", m_chiHcal));
  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
