/**
 *  @file   LCContent/src/LCTopologicalAssociation/BIBHitMergingAlgorithm.cc
 *
 *  @brief  Implementation of the BIB hit merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LCTopologicalAssociation/BIBHitMergingAlgorithm.h"

#include "LCUtility/KDTreeLinkerAlgoT.h"

#include <algorithm>
#include <array>
#include <limits>
#include <map>

using namespace pandora;

namespace lc_content {

BIBHitMergingAlgorithm::BIBHitMergingAlgorithm()
    : m_shouldUseCurrentClusterList(true), m_minHostClusterEnergy(2.f), m_maxRecombinationDistance(50.f),
      m_minClusterNeighbours(2), m_neighbourNLayers(2), m_neighbourCutDistanceSquared(25.f * 25.f),
      m_neighbourMaxSeparationSquared(1000.f * 1000.f), m_searchSafetyFactor(2.5f), m_layerSpanTolerance(2),
      m_timingCutEnabled(false), m_maxTimeDifference(1.f) {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BIBHitMergingAlgorithm::Run() {
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
      std::cout << "BIBHitMergingAlgorithm: Failed to obtain cluster list " << *iter << std::endl;
    }
  }

  // Eligible hosts: track-associated clusters, or clusters above the energy threshold.
  // Per host: pseudolayer span, energy-weighted mean hit time (for the optional timing
  // gate), and the kd-tree nodes of its non-isolated member hits.
  typedef std::pair<const CaloHit*, unsigned int> HitAndHost;
  std::vector<const Cluster*> hosts;
  std::vector<unsigned int> hostInnerLayer, hostOuterLayer;
  FloatVector hostMeanTime;
  std::vector<HitKDNode3D> memberNodes;
  std::array<float, 3> minpos{{0.f, 0.f, 0.f}}, maxpos{{0.f, 0.f, 0.f}};

  for (ClusterList::const_iterator iterJ = clusterList.begin(), iterJEnd = clusterList.end(); iterJ != iterJEnd;
       ++iterJ) {
    const Cluster* const pCluster = *iterJ;

    if (NULL == pCluster)
      continue;

    if (pCluster->GetAssociatedTrackList().empty() && (pCluster->GetHadronicEnergy() < m_minHostClusterEnergy))
      continue;

    const unsigned int hostIndex(hosts.size());
    hosts.push_back(pCluster);
    hostInnerLayer.push_back(pCluster->GetInnerPseudoLayer());
    hostOuterLayer.push_back(pCluster->GetOuterPseudoLayer());

    float timeEnergySum(0.f), energySum(0.f);
    const OrderedCaloHitList& orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end();
         iter != iterEnd; ++iter) {
      for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end();
           hitIter != hitIterEnd; ++hitIter) {
        const CaloHit* const pMemberHit = *hitIter;
        const CartesianVector& position(pMemberHit->GetPositionVector());
        memberNodes.emplace_back(HitAndHost(pMemberHit, hostIndex), position.GetX(), position.GetY(), position.GetZ());

        timeEnergySum += pMemberHit->GetTime() * pMemberHit->GetElectromagneticEnergy();
        energySum += pMemberHit->GetElectromagneticEnergy();

        if (1 == memberNodes.size()) {
          minpos[0] = position.GetX();
          minpos[1] = position.GetY();
          minpos[2] = position.GetZ();
          maxpos[0] = position.GetX();
          maxpos[1] = position.GetY();
          maxpos[2] = position.GetZ();
        } else {
          minpos[0] = std::min(position.GetX(), minpos[0]);
          minpos[1] = std::min(position.GetY(), minpos[1]);
          minpos[2] = std::min(position.GetZ(), minpos[2]);
          maxpos[0] = std::max(position.GetX(), maxpos[0]);
          maxpos[1] = std::max(position.GetY(), maxpos[1]);
          maxpos[2] = std::max(position.GetZ(), maxpos[2]);
        }
      }
    }

    hostMeanTime.push_back((energySum > 0.f) ? timeEnergySum / energySum : 0.f);
  }

  if (memberNodes.empty())
    return STATUS_CODE_SUCCESS;

  HitKDTree3D membersKdTree;
  KDTreeCube membersBoundingRegion(minpos[0], maxpos[0], minpos[1], maxpos[1], minpos[2], maxpos[2]);
  membersKdTree.build(memberNodes, membersBoundingRegion);

  // Loop over the available flagged hits and attach each to the host cluster it is
  // locally continuous with, if any. Attached hits become isolated hits of the host,
  // not members, so they never seed further attachments.
  const CaloHitList* pCaloHitList = NULL;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

  const float searchDistance(m_searchSafetyFactor * std::sqrt(m_neighbourCutDistanceSquared));
  std::vector<HitKDNode3D> found;

  for (CaloHitList::const_iterator hitIterI = pCaloHitList->begin(); hitIterI != pCaloHitList->end(); ++hitIterI) {
    const CaloHit* const pCaloHit = *hitIterI;

    if (!pCaloHit->IsPossibleBIB() || !PandoraContentApi::IsAvailable(*this, pCaloHit))
      continue;

    const CartesianVector& positionVector(pCaloHit->GetPositionVector());
    const float positionMagnitudeSquared(positionVector.GetMagnitudeSquared());
    const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());

    found.clear();
    KDTreeCube searchRegion = build_3d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance);
    membersKdTree.search(searchRegion, found);

    if (found.empty())
      continue;

    // Count cylinder neighbours and track the nearest member hit, per candidate host
    std::map<unsigned int, unsigned int> neighbourCounts;
    std::map<unsigned int, float> minDistanceSquareds;

    for (std::vector<HitKDNode3D>::const_iterator iterF = found.begin(), iterFEnd = found.end(); iterF != iterFEnd;
         ++iterF) {
      const CaloHit* const pMemberHit = iterF->data.first;
      const unsigned int hostIndex(iterF->data.second);

      const unsigned int memberLayer(pMemberHit->GetPseudoLayer());
      const unsigned int layerDifference((memberLayer > pseudoLayer) ? memberLayer - pseudoLayer
                                                                     : pseudoLayer - memberLayer);

      if (layerDifference > m_neighbourNLayers)
        continue;

      const CartesianVector positionDifference(positionVector - pMemberHit->GetPositionVector());
      const float separationSquared(positionDifference.GetMagnitudeSquared());

      if (separationSquared > m_neighbourMaxSeparationSquared)
        continue;

      std::map<unsigned int, float>::iterator distIter(minDistanceSquareds.find(hostIndex));

      if ((minDistanceSquareds.end() == distIter) || (separationSquared < distIter->second))
        minDistanceSquareds[hostIndex] = separationSquared;

      const CartesianVector crossProduct(positionVector.GetCrossProduct(positionDifference));

      if ((crossProduct.GetMagnitudeSquared() / positionMagnitudeSquared) < m_neighbourCutDistanceSquared)
        ++neighbourCounts[hostIndex];
    }

    // Best host: most cylinder neighbours, ties broken by nearest member hit
    const Cluster* pBestHostCluster(NULL);
    unsigned int bestNeighbourCount(0);
    float bestDistanceSquared(std::numeric_limits<float>::max());

    for (std::map<unsigned int, unsigned int>::const_iterator iterC = neighbourCounts.begin(),
                                                              iterCEnd = neighbourCounts.end();
         iterC != iterCEnd; ++iterC) {
      const unsigned int hostIndex(iterC->first);
      const unsigned int neighbourCount(iterC->second);

      if (neighbourCount < m_minClusterNeighbours)
        continue;

      const float distanceSquared(minDistanceSquareds[hostIndex]);

      if (distanceSquared > m_maxRecombinationDistance * m_maxRecombinationDistance)
        continue;

      if ((pseudoLayer + m_layerSpanTolerance < hostInnerLayer[hostIndex]) ||
          (pseudoLayer > hostOuterLayer[hostIndex] + m_layerSpanTolerance))
        continue;

      if (m_timingCutEnabled && (std::fabs(pCaloHit->GetTime() - hostMeanTime[hostIndex]) > m_maxTimeDifference))
        continue;

      if ((neighbourCount > bestNeighbourCount) ||
          ((neighbourCount == bestNeighbourCount) && (distanceSquared < bestDistanceSquared))) {
        pBestHostCluster = hosts[hostIndex];
        bestNeighbourCount = neighbourCount;
        bestDistanceSquared = distanceSquared;
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

StatusCode BIBHitMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {
  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "ShouldUseCurrentClusterList", m_shouldUseCurrentClusterList));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadVectorOfValues(xmlHandle, "AdditionalClusterListNames", m_additionalClusterListNames));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinHostClusterEnergy", m_minHostClusterEnergy));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
      XmlHelper::ReadValue(xmlHandle, "MaxRecombinationDistance", m_maxRecombinationDistance));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MinClusterNeighbours", m_minClusterNeighbours));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "NeighbourNLayers", m_neighbourNLayers));

  float neighbourCutDistance(std::sqrt(m_neighbourCutDistanceSquared));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "NeighbourCutDistance", neighbourCutDistance));
  m_neighbourCutDistanceSquared = neighbourCutDistance * neighbourCutDistance;

  float neighbourMaxSeparation(std::sqrt(m_neighbourMaxSeparationSquared));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "NeighbourMaxSeparation", neighbourMaxSeparation));
  m_neighbourMaxSeparationSquared = neighbourMaxSeparation * neighbourMaxSeparation;

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "SearchSafetyFactor", m_searchSafetyFactor));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "LayerSpanTolerance", m_layerSpanTolerance));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "TimingCutEnabled", m_timingCutEnabled));

  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                  XmlHelper::ReadValue(xmlHandle, "MaxTimeDifference", m_maxTimeDifference));

  return STATUS_CODE_SUCCESS;
}

} // namespace lc_content
