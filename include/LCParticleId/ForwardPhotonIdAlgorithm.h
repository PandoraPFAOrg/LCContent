/**
 *  @file   LCContent/include/LCParticleId/ForwardPhotonIdAlgorithm.h
 *
 *  @brief  Cut-based forward-photon tagging for HCAL-only clusters.
 *
 *  A trackless HCAL-only cluster is tagged as a photon (ParticleId = PHOTON) iff
 *      |eta|        > MinAbsEta         (forward region)      [evaluated first]
 *      C / S        > MinCS             (Cherenkov-rich, EM-like)
 *      coreFraction > MinCoreFraction   (laterally compact)
 *  where C, S are the Cherenkov/scintillation HCAL hit-energy sums, eta is the
 *  pseudorapidity of the W0 log-weighted centroid, and coreFraction is the
 *  dual-readout energy within a CoreMm cone (half-angle atan(CoreMm/r_centroid))
 *  of the centroid divided by the cluster's total dual-readout energy.
 *
 *  The |eta| cut is applied first: it is cheap and rejects the bulk of (central)
 *  clusters before the C/S and core-fraction passes. Clusters that fail any cut
 *  are left untouched. The tag decouples forward HCAL photons before the
 *  satellite / neutral-primary steps downstream.
 */
#ifndef LC_FORWARD_PHOTON_ID_ALGORITHM_H
#define LC_FORWARD_PHOTON_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"

#include <vector>

namespace lc_content {

class ForwardPhotonIdAlgorithm : public pandora::Algorithm {
public:
  class Factory : public pandora::AlgorithmFactory {
  public:
    pandora::Algorithm *CreateAlgorithm() const override {
      return new ForwardPhotonIdAlgorithm();
    }
  };

  ForwardPhotonIdAlgorithm() = default;

  /// A single HCAL hit: position, hadronic-scale energy and Cherenkov flag.
  /// Public so it can name the helper-method signatures defined out-of-line.
  struct HcalHit { float x, y, z, e; bool isCher; };

private:
  pandora::StatusCode Run() override;
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

  /// Collect the cluster's HCAL hits and their total energy. Returns false if
  /// the cluster carries any ECAL hit (not HCAL-only) or has no HCAL energy.
  bool CollectHcalHits(const pandora::Cluster *const pCluster,
                       std::vector<HcalHit> &hits, float &eTotal) const;

  /// W0 log-weighted centroid of the HCAL hits. Returns false if degenerate.
  bool W0Centroid(const std::vector<HcalHit> &hits, const float eTotal,
                  pandora::CartesianVector &centroid) const;

  /// Cherenkov-to-scintillation energy ratio C/S; returns a negative sentinel
  /// when S <= 0 (pure-Cherenkov / empty), which the caller treats as a fail.
  float CherToScintRatio(const std::vector<HcalHit> &hits) const;

  /// DR core-energy fraction within a CoreMm cone of the centroid direction.
  float CoreFraction(const std::vector<HcalHit> &hits,
                     const pandora::CartesianVector &centroid,
                     const float rCentroid) const;

  /// Dual-readout HCAL energy of a hit set: (S - ChiHcal * C) / (1 - ChiHcal).
  float DRcorrHcal(const std::vector<HcalHit> &hits) const;

  float m_minCS       = 0.5f;    ///< min Cherenkov/scintillation ratio C/S
  float m_minAbsEta   = 2.4f;    ///< min |eta| of the cluster centroid
  float m_coreMm      = 25.f;    ///< core cone lateral size [mm]
  float m_minCoreFrac = 0.6f;    ///< min DR core-energy fraction
  float m_w0          = 4.6f;    ///< log-weight offset for the centroid
  float m_chiHcal     = 0.31f;   ///< HCAL dual-readout chi factor
};

} // namespace lc_content

#endif // LC_FORWARD_PHOTON_ID_ALGORITHM_H
