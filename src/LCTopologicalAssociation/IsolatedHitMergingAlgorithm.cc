/**
 *  @file   LCContent/src/LCTopologicalAssociation/IsolatedHitMergingAlgorithm.cc
 *
 *  @brief  Implementation of the isolated hit merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include "LCHelpers/SortingHelper.h"

#include "LCTopologicalAssociation/IsolatedHitMergingAlgorithm.h"

#include <algorithm>
#include <cmath>

using namespace pandora;

namespace lc_content {

IsolatedHitMergingAlgorithm::IsolatedHitMergingAlgorithm()
    : m_shouldUseCurrentClusterList(true), m_minHitsInCluster(4), m_maxRecombinationDistance(250.f),
      m_minCosOpeningAngle(0.f), m_useCorrectedHadronicEnergy(false), m_ignorePhotons(false), m_ignoreCharged(false),
      m_useAngularDistance(false), m_shouldExcludeBIBHits(false) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::Run() {
  // HACK
  const ClusterList* pInputClusterList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputClusterList));

  // Read specified lists of input clusters
  ClusterList clusterList;

  if (m_shouldUseCurrentClusterList) {
    const ClusterList* pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
  }

  for (StringVector::const_iterator iter = m_additionalClusterListNames.begin(),
                                    iterEnd = m_additionalClusterListNames.end();
       iter != iterEnd; ++iter) {
    const ClusterList* pClusterList = NULL;

    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, *iter, pClusterList)) {
      clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    } else {
      std::cout << "IsolatedHitMergingAlgorithm: Failed to obtain cluster list " << *iter << std::endl;
    }
  }

  // Create a vector of input clusters, ordered by inner layer
  ClusterVector clusterVector(clusterList.begin(), clusterList.end());
  std::sort(clusterVector.begin(), clusterVector.end(), SortingHelper::SortClustersByInnerLayer);

  // Build cache for remaining clusters
  std::vector<ClusterCache> clusterCaches;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildClusterCache(clusterVector, clusterCaches));

  // FIRST PART - find "small" clusters, below threshold number of calo hits, delete them and associate hits with other
  // clusters
  for (std::size_t iClus = 0; iClus < clusterVector.size(); ++iClus) {
    const Cluster* const pClusterToDelete = clusterVector[iClus];

    if (NULL == pClusterToDelete)
      continue;

    const unsigned int nCaloHits(pClusterToDelete->GetNCaloHits());

    if (nCaloHits > m_minHitsInCluster)
      continue;

    if (pInputClusterList->end() == std::find(pInputClusterList->begin(), pInputClusterList->end(), pClusterToDelete))
      continue;

    CaloHitList caloHitList;
    pClusterToDelete->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pClusterToDelete));
    // Null both the vector entry and the parallel cache entry.
    clusterVector[iClus] = NULL;
    clusterCaches[iClus].pCluster = NULL;

    // Redistribute hits that used to be in cluster I amongst other clusters
    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end();
         hitIter != hitIterEnd; ++hitIter) {
      const CaloHit* const pCaloHit = *hitIter;

      if (m_shouldExcludeBIBHits && pCaloHit->IsPossibleBIB())
        continue;

      const Cluster* pBestHostCluster(NULL);
      float bestHostClusterEnergy(0.);
      float minDistance(m_maxRecombinationDistance);

      // Find the most appropriate cluster for this newly-available hit
      for (const ClusterCache& cache : clusterCaches) {
        if (!cache.pCluster)
          continue;

        if (cache.pCluster->GetNCaloHits() < nCaloHits)
          continue;

        if (m_ignoreCharged && cache.isCharged)
          continue;

        float distance = GetDistanceToHit(cache, pCaloHit);

        // In event of equidistant host candidates, choose highest energy cluster
        if ((distance < minDistance) || (distance == minDistance && cache.energy > bestHostClusterEnergy)) {
          minDistance = distance;
          pBestHostCluster = cache.pCluster;
          bestHostClusterEnergy = cache.energy;
        }
      }

      if (NULL != pBestHostCluster) {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                 PandoraContentApi::AddIsolatedToCluster(*this, pBestHostCluster, pCaloHit));
      }
    }
  }

  // SECOND PART - loop over isolated hits and associate each with the best cluster.
  const CaloHitList* pCaloHitList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

  for (const auto* pCaloHit : *pCaloHitList) {
    if (!pCaloHit->IsIsolated() || !PandoraContentApi::IsAvailable(*this, pCaloHit))
      continue;

    if (m_shouldExcludeBIBHits && pCaloHit->IsPossibleBIB())
      continue;

    const Cluster* pBestHostCluster(NULL);
    float bestHostClusterEnergy(0.f);
    float minDistance(m_maxRecombinationDistance);

    for (const ClusterCache& cache : clusterCaches) {
      if (!cache.pCluster)
        continue;
      if (m_ignorePhotons && cache.isEm && !cache.isCharged)
        continue;
      if (m_ignoreCharged && cache.isCharged)
        continue;

      float distance = GetDistanceToHit(cache, pCaloHit);

      if ((distance < minDistance) || (distance == minDistance && cache.energy > bestHostClusterEnergy)) {
        minDistance = distance;
        pBestHostCluster = cache.pCluster;
        bestHostClusterEnergy = cache.energy;
      }
    }

    if (NULL != pBestHostCluster) {
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                               PandoraContentApi::AddIsolatedToCluster(*this, pBestHostCluster, pCaloHit));
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::BuildClusterCache(const ClusterVector& clusterVector,
                                                          std::vector<ClusterCache>& clusterCaches) const {
  clusterCaches.clear();
  clusterCaches.reserve(clusterVector.size());

  for (const Cluster* const pCluster : clusterVector) {
    ClusterCache cache;
    cache.pCluster = pCluster; // may be nullptr if already deleted in Part 1

    if (pCluster) {
      cache.centroid = pCluster->GetCentroid(pCluster->GetInnerPseudoLayer());
      cache.direction = pCluster->GetInitialDirection();
      cache.isCharged = !pCluster->GetAssociatedTrackList().empty();
      cache.isEm = pCluster->GetParticleId() == pandora::PHOTON;
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusterEnergy(pCluster, cache.energy));
    } else {
      cache.energy = 0.f;
      cache.isEm = false;
      cache.isCharged = false;
    }

    clusterCaches.push_back(cache);
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float IsolatedHitMergingAlgorithm::GetDistanceToHit(const ClusterCache& cache, const CaloHit* const pCaloHit) const {
  // Apply simple preselection using cosine of opening angle between the hit and cluster directions
  if (pCaloHit->GetExpectedDirection().GetCosOpeningAngle(cache.direction) < m_minCosOpeningAngle) {
    return std::numeric_limits<float>::max();
  }

  const CartesianVector& hitPosition(pCaloHit->GetPositionVector());
  float distance = std::numeric_limits<float>::max();

  // Angular mode: use (1 - cosA) as the distance metric.
  if (m_useAngularDistance) {
    float cosA = cache.centroid.GetCosOpeningAngle(hitPosition);
    cosA = std::max(-1.f, std::min(1.f, cosA));
    distance = 1.f - cosA;
  } else {
    const OrderedCaloHitList& orderedCaloHitList(cache.pCluster->GetOrderedCaloHitList());
    float minDistanceSquared = std::numeric_limits<float>::max();

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end();
         iter != iterEnd; ++iter) {
      const float distanceSquared = (cache.pCluster->GetCentroid(iter->first) - hitPosition).GetMagnitudeSquared();

      if (distanceSquared < minDistanceSquared)
        minDistanceSquared = distanceSquared;
    }

    distance = std::sqrt(minDistanceSquared);
  }

  return distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::GetClusterEnergy(const Cluster* const pCluster, float& energy) const {
  if (!m_useCorrectedHadronicEnergy) {
    energy = pCluster->GetHadronicEnergy();
    return STATUS_CODE_SUCCESS;
  }

  const pandora::EnergyCorrections* const pEnergyCorrections(
      PandoraContentApi::GetPlugins(*this)->GetEnergyCorrections());
  if (!pEnergyCorrections)
    return STATUS_CODE_FAILURE;

  float corrEm(0.f), corrHad(0.f);
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                           pEnergyCorrections->MakeEnergyCorrections(pCluster, corrEm, corrHad));
  energy = corrHad;

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IsolatedHitMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ShouldUseCurrentClusterList", m_shouldUseCurrentClusterList));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "ShouldExcludeBIBHits", m_shouldExcludeBIBHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadVectorOfValues(xmlHandle, "AdditionalClusterListNames", m_additionalClusterListNames));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinHitsInCluster", m_minHitsInCluster));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "MaxRecombinationDistance", m_maxRecombinationDistance));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinCosOpeningAngle", m_minCosOpeningAngle));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "UseCorrectedHadronicEnergy", m_useCorrectedHadronicEnergy));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "UseAngularDistance", m_useAngularDistance));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "IgnorePhotons", m_ignorePhotons));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "IgnoreCharged", m_ignoreCharged));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
