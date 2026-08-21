/**
 *  @file   LCContent/src/LCHelpers/ClusterProximityHelper.cc
 *
 *  @brief  Implementation of the cluster proximity helper classes.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/ClusterProximityHelper.h"

#include <algorithm>
#include <utility>

using namespace pandora;

namespace lc_content {

ClusterBoundingBox::ClusterBoundingBox() : m_min(0.f, 0.f, 0.f), m_max(0.f, 0.f, 0.f), m_isInitialized(false) {}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterBoundingBox::ClusterBoundingBox(const Cluster* const pCluster)
    : m_min(0.f, 0.f, 0.f), m_max(0.f, 0.f, 0.f), m_isInitialized(false) {
  for (const auto& layerEntry : pCluster->GetOrderedCaloHitList())
    this->Enclose(*(layerEntry.second));
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterBoundingBox::ClusterBoundingBox(const CaloHitList& caloHitList)
    : m_min(0.f, 0.f, 0.f), m_max(0.f, 0.f, 0.f), m_isInitialized(false) {
  this->Enclose(caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterBoundingBox::Enclose(const CaloHitList& caloHitList) {
  for (const CaloHit* const pCaloHit : caloHitList) {
    const CartesianVector& position(pCaloHit->GetPositionVector());

    if (!m_isInitialized) {
      m_min = position;
      m_max = position;
      m_isInitialized = true;
    } else {
      m_min.SetValues(std::min(m_min.GetX(), position.GetX()), std::min(m_min.GetY(), position.GetY()),
                      std::min(m_min.GetZ(), position.GetZ()));
      m_max.SetValues(std::max(m_max.GetX(), position.GetX()), std::max(m_max.GetY(), position.GetY()),
                      std::max(m_max.GetZ(), position.GetZ()));
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterBoundingBox& ClusterContactCache::GetBoundingBox(const Cluster* const pCluster) {
  const auto iter = m_boundingBoxes.find(pCluster);

  if (m_boundingBoxes.end() != iter)
    return iter->second;

  return m_boundingBoxes.emplace(pCluster, ClusterBoundingBox(pCluster)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterLayerBoundingBoxVector& ClusterContactCache::GetLayerBoundingBoxes(const Cluster* const pCluster) {
  const auto iter = m_layerBoundingBoxes.find(pCluster);

  if (m_layerBoundingBoxes.end() != iter)
    return iter->second;

  const OrderedCaloHitList& orderedCaloHitList(pCluster->GetOrderedCaloHitList());

  ClusterLayerBoundingBoxVector layerBoundingBoxes;
  layerBoundingBoxes.reserve(orderedCaloHitList.size());

  for (const auto& layerEntry : orderedCaloHitList)
    layerBoundingBoxes.emplace_back(*(layerEntry.second));

  return m_layerBoundingBoxes.emplace(pCluster, std::move(layerBoundingBoxes)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterContactCache::RecordMerge(const Cluster* const pParentCluster, const Cluster* const pDaughterCluster) {
  // The parent's hits have changed and the daughter is gone; every other box still holds.
  m_boundingBoxes.erase(pParentCluster);
  m_boundingBoxes.erase(pDaughterCluster);
  m_layerBoundingBoxes.erase(pParentCluster);
  m_layerBoundingBoxes.erase(pDaughterCluster);
  m_watermarks.erase(pDaughterCluster);

  m_merges.emplace_back(pParentCluster, pDaughterCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterVector ClusterContactCache::GetClustersChangedSince(const Cluster* const pCluster) const {
  const auto iter = m_watermarks.find(pCluster);
  const unsigned int watermark((m_watermarks.end() != iter) ? iter->second : 0);

  ClusterVector changedClusters;

  // A vector, and not the ClusterSet below, because ClusterSet is unordered and the caller needs a
  // reproducible order. The set answers "seen already" only - a daughter left alone for many passes has
  // many merges to catch up on, and a linear scan per entry would make that quadratic.
  ClusterSet seenClusters;

  for (unsigned int iMerge = watermark, nMerges = this->GetNMerges(); iMerge < nMerges; ++iMerge) {
    const ClusterMerge& merge(m_merges.at(iMerge));

    for (const Cluster* const pChangedCluster : {merge.first, merge.second}) {
      if (seenClusters.insert(pChangedCluster).second)
        changedClusters.push_back(pChangedCluster);
    }
  }

  return changedClusters;
}

} // namespace lc_content
