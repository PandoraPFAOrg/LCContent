/**
 *  @file   LCContent/src/LCHelpers/SortingHelper.cc
 *
 *  @brief  Implementation of the sorting helper class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/StatusCodes.h"

#include "LCHelpers/SortingHelper.h"

#include <cmath>
#include <limits>

using namespace pandora;

namespace lc_content {

bool SortingHelper::SortClustersByNHits(const Cluster* const pLhs, const Cluster* const pRhs) {
  // NHits
  const unsigned int nCaloHitsLhs(pLhs->GetNCaloHits()), nCaloHitsRhs(pRhs->GetNCaloHits());

  if (nCaloHitsLhs != nCaloHitsRhs)
    return (nCaloHitsLhs > nCaloHitsRhs);

  // Track seeds
  if ((0 == nCaloHitsLhs) && (0 == nCaloHitsRhs)) {
    const float trackEnergyLhs(pLhs->IsTrackSeeded() ? pLhs->GetTrackSeed()->GetEnergyAtDca() : 0.f);
    const float trackEnergyRhs(pRhs->IsTrackSeeded() ? pRhs->GetTrackSeed()->GetEnergyAtDca() : 0.f);

    if (std::fabs(trackEnergyLhs - trackEnergyRhs) > std::numeric_limits<float>::epsilon())
      return (trackEnergyLhs > trackEnergyRhs);
  }

  // Energy
  const float energyLhs(pLhs->GetHadronicEnergy()), energyRhs(pRhs->GetHadronicEnergy());

  if (std::fabs(energyLhs - energyRhs) > std::numeric_limits<float>::epsilon())
    return (energyLhs > energyRhs);

  // Energy in isolated hits
  const float isolatedEnergyLhs(pLhs->GetIsolatedHadronicEnergy()),
      isolatedEnergyRhs(pRhs->GetIsolatedHadronicEnergy());

  if (std::fabs(isolatedEnergyLhs - isolatedEnergyRhs) > std::numeric_limits<float>::epsilon())
    return (isolatedEnergyLhs > isolatedEnergyRhs);

  // Final attempt to distinguish
  if ((nCaloHitsLhs > 0) && (nCaloHitsRhs > 0)) {
    const CaloHit* const pFirstHitLhs((pLhs->GetOrderedCaloHitList().begin())->second->front());
    const CaloHit* const pFirstHitRhs((pRhs->GetOrderedCaloHitList().begin())->second->front());

    return (*pFirstHitLhs < *pFirstHitRhs);
  }

  throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SortingHelper::SortClustersByInnerLayer(const Cluster* const pLhs, const Cluster* const pRhs) {
  const unsigned int nCaloHitsLhs(pLhs->GetNCaloHits()), nCaloHitsRhs(pRhs->GetNCaloHits());

  if ((nCaloHitsLhs > 0) && (nCaloHitsRhs > 0)) {
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer()), innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if (innerLayerLhs != innerLayerRhs)
      return (innerLayerLhs < innerLayerRhs);
  }

  return SortingHelper::SortClustersByNHits(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SortingHelper::SortPfosByEnergy(const ParticleFlowObject* const pLhs, const ParticleFlowObject* const pRhs) {
  return (pLhs->GetEnergy() > pRhs->GetEnergy());
}

} // namespace lc_content
