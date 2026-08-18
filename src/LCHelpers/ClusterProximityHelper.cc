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

namespace lc_content
{

ClusterBoundingBox::ClusterBoundingBox() :
    m_min{0.f, 0.f, 0.f},
    m_max{0.f, 0.f, 0.f},
    m_isInitialized(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterBoundingBox::ClusterBoundingBox(const Cluster *const pCluster) :
    m_min{0.f, 0.f, 0.f},
    m_max{0.f, 0.f, 0.f},
    m_isInitialized(false)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
        this->Enclose(*(iter->second));
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterBoundingBox::ClusterBoundingBox(const CaloHitList &caloHitList) :
    m_min{0.f, 0.f, 0.f},
    m_max{0.f, 0.f, 0.f},
    m_isInitialized(false)
{
    this->Enclose(caloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterBoundingBox::Enclose(const CaloHitList &caloHitList)
{
    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CartesianVector &position((*hitIter)->GetPositionVector());
        const float coordinates[3] = {position.GetX(), position.GetY(), position.GetZ()};

        if (!m_isInitialized)
        {
            for (unsigned int i = 0; i < 3; ++i)
                m_min[i] = m_max[i] = coordinates[i];

            m_isInitialized = true;
        }
        else
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
                m_min[i] = std::min(m_min[i], coordinates[i]);
                m_max[i] = std::max(m_max[i], coordinates[i]);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterBoundingBox &ClusterContactCache::GetBoundingBox(const Cluster *const pCluster)
{
    ClusterToBoundingBoxMap::const_iterator iter = m_boundingBoxes.find(pCluster);

    if (m_boundingBoxes.end() != iter)
        return iter->second;

    return m_boundingBoxes.emplace(pCluster, ClusterBoundingBox(pCluster)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ClusterLayerBoundingBoxVector &ClusterContactCache::GetLayerBoundingBoxes(const Cluster *const pCluster)
{
    ClusterToLayerBoundingBoxMap::const_iterator iter = m_layerBoundingBoxes.find(pCluster);

    if (m_layerBoundingBoxes.end() != iter)
        return iter->second;

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    ClusterLayerBoundingBoxVector layerBoundingBoxes;
    layerBoundingBoxes.reserve(orderedCaloHitList.size());

    for (OrderedCaloHitList::const_iterator hitIter = orderedCaloHitList.begin(), hitIterEnd = orderedCaloHitList.end(); hitIter != hitIterEnd; ++hitIter)
        layerBoundingBoxes.emplace_back(*(hitIter->second));

    return m_layerBoundingBoxes.emplace(pCluster, std::move(layerBoundingBoxes)).first->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterContactCache::RecordMerge(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster)
{
    // The parent's hits have changed and the daughter is gone; every other box still holds.
    m_boundingBoxes.erase(pParentCluster);
    m_boundingBoxes.erase(pDaughterCluster);
    m_layerBoundingBoxes.erase(pParentCluster);
    m_layerBoundingBoxes.erase(pDaughterCluster);
    m_watermarks.erase(pDaughterCluster);

    m_merges.emplace_back(pParentCluster, pDaughterCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterContactCache::GetClustersChangedSince(const Cluster *const pCluster, ClusterVector &changedClusters) const
{
    ClusterToWatermarkMap::const_iterator iter = m_watermarks.find(pCluster);
    const unsigned int watermark((m_watermarks.end() != iter) ? iter->second : 0);

    // Deterministic order, from a set used only to answer "seen already" - a daughter left alone for many
    // passes has many merges to catch up on, and a linear scan per entry makes that quadratic.
    ClusterSet seenClusters;

    for (unsigned int iMerge = watermark, nMerges = this->GetNMerges(); iMerge < nMerges; ++iMerge)
    {
        const ClusterMerge &merge(m_merges.at(iMerge));

        for (const Cluster *const pChangedCluster : {merge.first, merge.second})
        {
            if (seenClusters.insert(pChangedCluster).second)
                changedClusters.push_back(pChangedCluster);
        }
    }
}

} // namespace lc_content
