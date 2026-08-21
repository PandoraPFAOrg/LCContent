/**
 *  @file   LCContent/include/LCFragmentRemoval/FragmentRemovalBaseAlgorithm.h
 *
 *  @brief  Header file for the fragment removal base algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_FRAGMENT_REMOVAL_BASE_ALGORITHM_H
#define LC_FRAGMENT_REMOVAL_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "LCHelpers/ClusterProximityHelper.h"
#include "LCHelpers/FragmentRemovalHelper.h"

#include <algorithm>
#include <limits>
#include <map>
#include <vector>

namespace lc_content {

/**
 *  @brief  FragmentRemovalBaseAlgorithm class, the machinery shared by the fragment removal algorithms
 *
 *          Each of those algorithms merges one cluster pair per pass, then reconsiders the clusters the
 *          merge disturbed and goes round again. What differs between them is which clusters may be a
 *          daughter, which may be a parent, what makes a contact worth storing and how the winning pair
 *          is chosen; the pass loop and the bookkeeping around it are the same in all of them, and live
 *          here. The differences are supplied by the hooks below.
 *
 *          The template parameter is the cluster contact type, which fixes the contact parameters the
 *          algorithm reads and the shape of the map it fills.
 */
template <typename CONTACT>
class FragmentRemovalBaseAlgorithm : public pandora::Algorithm {
public:
  using ContactParameters = typename CONTACT::Parameters;
  using ContactVector = std::vector<CONTACT>;
  using ContactMap = std::map<const pandora::Cluster*, ContactVector>;

  /**
   *  @brief  Default constructor
   */
  FragmentRemovalBaseAlgorithm();

protected:
  pandora::StatusCode Run();

  /**
   *  @brief  Read the settings shared by every fragment removal algorithm
   *
   *          A derived class calls this from its own ReadSettings, then reads whatever else it needs.
   *          Defaults belong in the derived constructor, since they differ from algorithm to algorithm.
   *
   *  @param  xmlHandle the relevant xml handle
   */
  pandora::StatusCode ReadCommonSettings(const pandora::TiXmlHandle xmlHandle);

  /**
   *  @brief  Get the clusters this algorithm considers, from the current cluster list
   *
   *          The order of the returned list is the order contact vectors are built in, and the merging
   *          candidate search breaks exact ties on the first entry it meets, so it has to be stable.
   *          The default keeps the current list as it stands.
   *
   *  @param  inputClusterList the current cluster list
   *  @param  clusterList to receive the clusters to consider
   */
  virtual void GetRelevantClusterList(const pandora::ClusterList& inputClusterList,
                                      pandora::ClusterList& clusterList) const;

  /**
   *  @brief  Whether a cluster may be a daughter candidate, beyond the hit count and hadronic energy cuts
   *
   *  @param  pDaughterCluster address of the daughter candidate cluster
   *
   *  @return boolean
   */
  virtual bool IsCandidateDaughter(const pandora::Cluster* const pDaughterCluster) const = 0;

  /**
   *  @brief  Whether a cluster may be a parent candidate for the provided daughter
   *
   *          Called for every ordered pair of distinct clusters, so it must not depend on anything but
   *          its two clusters: the incremental update recomputes a contact only when one of the two has
   *          been changed by a merge.
   *
   *  @param  pDaughterCluster address of the daughter candidate cluster
   *  @param  pParentCluster address of the parent candidate cluster
   *
   *  @return boolean
   */
  virtual bool IsCandidateParent(const pandora::Cluster* const pDaughterCluster,
                                 const pandora::Cluster* const pParentCluster) const = 0;

  /**
   *  @brief  Whether candidate parent and daughter clusters are sufficiently in contact to warrant further
   * investigation
   *
   *          Only reached for a contact whose closest hit-hit separation is already within
   *          m_contactCutMaxDistance: PassesContactCuts applies that cut before dispatching here, so an
   *          implementation neither has to repeat it nor is able to relax it. See PassesContactCuts for
   *          why it is not an implementation's to make.
   *
   *  @param  contact the cluster contact
   *
   *  @return boolean
   */
  virtual bool PassesClusterContactCuts(const CONTACT& contact) const = 0;

  /**
   *  @brief  Find the best candidate parent and daughter clusters for fragment removal merging
   *
   *  @param  contactMap the populated cluster contact map
   *  @param  pBestParentCluster to receive the address of the best parent cluster candidate
   *  @param  pBestDaughterCluster to receive the address of the best daughter cluster candidate
   */
  virtual pandora::StatusCode GetClusterMergingCandidates(const ContactMap& contactMap,
                                                          const pandora::Cluster*& pBestParentCluster,
                                                          const pandora::Cluster*& pBestDaughterCluster) = 0;

  /**
   *  @brief  Act on the parent cluster once a merge has been applied
   *
   *          This is the one hook allowed to change a cluster. It is safe for it to do so because the
   *          merge parent is already in the set of clusters the merge changed, so every contact involving
   *          it is recomputed on the next pass regardless. The default does nothing.
   *
   *  @param  pBestParentCluster address of the parent cluster of the merge just applied
   */
  virtual pandora::StatusCode PostMergeAction(const pandora::Cluster* const pBestParentCluster);

  /**
   *  @brief  Whether a contact is worth storing
   *
   *          The max distance cut lives here, and not in the PassesClusterContactCuts implementations,
   *          because CouldPassClusterContactCuts skips whole pairs on the strength of it. Were a
   *          subclass free to drop or relax it, that subclass would never see the pairs it had just
   *          decided to keep, and nothing would say so. Holding the cut in the one place every contact
   *          passes through is what makes the skip provably safe rather than safe by convention.
   *
   *  @param  contact the cluster contact
   *
   *  @return boolean
   */
  bool PassesContactCuts(const CONTACT& contact) const;

  /**
   *  @brief  Whether a cluster contact between the two clusters could possibly pass PassesContactCuts
   *
   *          Decided from cached cluster properties alone, so that pairs which cannot contribute never
   *          enter the hit loops inside the contact constructor. A false return is a guarantee, not a
   *          heuristic, and rests on two things, each of which is now a property of code rather than a
   *          request to whoever writes the next subclass:
   *             - PassesContactCuts discards a contact whose closest hit-hit separation exceeds
   *               m_contactCutMaxDistance before any subclass term is consulted, and no subclass term
   *               can put such a contact back.
   *             - a contact whose clusters fail ClusterContact::PassesDirectionPreselection reports
   *               that separation as float max, because the constructor then skips the hit loop
   *               altogether. This asks that same function, so the two cannot answer differently.
   *
   *  @param  pDaughterCluster address of the daughter candidate cluster
   *  @param  daughterBoundingBox bounding box of the daughter candidate cluster
   *  @param  pParentCluster address of the parent candidate cluster
   *  @param  contactCache the bounding boxes and merge history carried across passes
   *
   *  @return boolean
   */
  bool CouldPassClusterContactCuts(const pandora::Cluster* const pDaughterCluster,
                                   const ClusterBoundingBox& daughterBoundingBox,
                                   const pandora::Cluster* const pParentCluster,
                                   ClusterContactCache& contactCache) const;

  ContactParameters m_contactParameters; ///< The cluster contact parameters

  unsigned int m_nMaxPasses;          ///< Maximum number of passes over cluster contact information
  unsigned int m_minDaughterCaloHits; ///< Min number of calo hits in daughter candidate clusters
  float m_minDaughterHadronicEnergy;  ///< Min hadronic energy for daughter candidate clusters
  float m_contactCutMaxDistance;      ///< Max distance between closest hits to store cluster contact info

private:
  /**
   *  @brief  Get cluster contact map, linking each daughter candidate cluster to a list of parent candidates and
   * describing the proximity/contact between each pairing
   *
   *  @param  isFirstPass whether this is the first call to GetClusterContactMap
   *  @param  affectedClusters list of those clusters affected by previous cluster merging, for which contact details
   * must be updated
   *  @param  contactCache the bounding boxes and merge history carried across passes
   *  @param  contactMap to receive the populated cluster contact map
   */
  pandora::StatusCode GetClusterContactMap(bool& isFirstPass, const pandora::ClusterSet& affectedClusters,
                                           ClusterContactCache& contactCache, ContactMap& contactMap) const;

  /**
   *  @brief  Bring one daughter's contact vector up to date after a merge, without rebuilding it
   *
   *          A merge changes exactly two clusters, so all but at most two of a daughter's contacts survive it
   *          unchanged. Rebuilding the vector recomputes every one of them; this recomputes the ones that moved.
   *
   *  @param  pDaughterCluster address of the daughter candidate cluster
   *  @param  changedClusters the clusters changed by the merges this daughter's contacts have not yet seen
   *  @param  changedClusterSet the same clusters, for membership tests
   *  @param  clusterToIndex position of each cluster in the relevant cluster list
   *  @param  contactCache the bounding boxes and merge history carried across passes
   *  @param  contactMap the cluster contact map to update
   */
  void UpdateClusterContacts(const pandora::Cluster* const pDaughterCluster,
                             const pandora::ClusterVector& changedClusters,
                             const pandora::ClusterSet& changedClusterSet, const ClusterToIndexMap& clusterToIndex,
                             ClusterContactCache& contactCache, ContactMap& contactMap) const;

  /**
   *  @brief  Get the list of clusters for which cluster contact information will be affected by a specified cluster
   * merge
   *
   *  @param  contactMap the cluster contact map
   *  @param  pBestParentCluster address of the parent cluster to be merged
   *  @param  pBestDaughterCluster address of the daughter cluster to be merged
   *  @param  affectedClusters to receive the list of affected clusters
   */
  pandora::StatusCode GetAffectedClusters(const ContactMap& contactMap,
                                          const pandora::Cluster* const pBestParentCluster,
                                          const pandora::Cluster* const pBestDaughterCluster,
                                          pandora::ClusterSet& affectedClusters) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
inline FragmentRemovalBaseAlgorithm<CONTACT>::FragmentRemovalBaseAlgorithm()
    : m_contactParameters(), m_nMaxPasses(std::numeric_limits<unsigned int>::max()), m_minDaughterCaloHits(0),
      m_minDaughterHadronicEnergy(0.f), m_contactCutMaxDistance(0.f) {}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
inline void FragmentRemovalBaseAlgorithm<CONTACT>::GetRelevantClusterList(const pandora::ClusterList& inputClusterList,
                                                                          pandora::ClusterList& clusterList) const {
  clusterList = inputClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
inline pandora::StatusCode
FragmentRemovalBaseAlgorithm<CONTACT>::PostMergeAction(const pandora::Cluster* const /*pBestParentCluster*/) {
  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
inline bool FragmentRemovalBaseAlgorithm<CONTACT>::PassesContactCuts(const CONTACT& contact) const {
  if (contact.GetDistanceToClosestHit() > m_contactCutMaxDistance)
    return false;

  return this->PassesClusterContactCuts(contact);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
inline bool FragmentRemovalBaseAlgorithm<CONTACT>::CouldPassClusterContactCuts(
    const pandora::Cluster* const pDaughterCluster, const ClusterBoundingBox& daughterBoundingBox,
    const pandora::Cluster* const pParentCluster, ClusterContactCache& contactCache) const {
  // See the header: both tests are rejections PassesContactCuts would make too, taken before the hits are
  // touched rather than after. The first is asked of the function the contact constructor itself asks, and
  // the second of the same m_contactCutMaxDistance the cut above uses.
  if (!ClusterContact::PassesDirectionPreselection(pDaughterCluster, pParentCluster, m_contactParameters))
    return false;

  return !daughterBoundingBox.IsSeparatedFrom(contactCache.GetBoundingBox(pParentCluster), m_contactCutMaxDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
pandora::StatusCode FragmentRemovalBaseAlgorithm<CONTACT>::Run() {
  unsigned int nPasses(0);
  bool isFirstPass(true), shouldRecalculate(true);

  pandora::ClusterSet affectedClusters;
  ContactMap contactMap;

  // Bounding boxes and the record of which clusters each merge changed, both reused across passes.
  ClusterContactCache contactCache;

  while ((nPasses++ < m_nMaxPasses) && shouldRecalculate) {
    shouldRecalculate = false;
    const pandora::Cluster *pBestParentCluster(NULL), *pBestDaughterCluster(NULL);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             this->GetClusterContactMap(isFirstPass, affectedClusters, contactCache, contactMap));

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                             this->GetClusterMergingCandidates(contactMap, pBestParentCluster, pBestDaughterCluster));

    if ((NULL != pBestParentCluster) && (NULL != pBestDaughterCluster)) {
      PANDORA_RETURN_RESULT_IF(
          pandora::STATUS_CODE_SUCCESS, !=,
          this->GetAffectedClusters(contactMap, pBestParentCluster, pBestDaughterCluster, affectedClusters));

      contactMap.erase(contactMap.find(pBestDaughterCluster));
      shouldRecalculate = true;

      contactCache.RecordMerge(pBestParentCluster, pBestDaughterCluster);

      PANDORA_RETURN_RESULT_IF(
          pandora::STATUS_CODE_SUCCESS, !=,
          PandoraContentApi::MergeAndDeleteClusters(*this, pBestParentCluster, pBestDaughterCluster));

      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->PostMergeAction(pBestParentCluster));
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
pandora::StatusCode FragmentRemovalBaseAlgorithm<CONTACT>::GetClusterContactMap(
    bool& isFirstPass, const pandora::ClusterSet& affectedClusters, ClusterContactCache& contactCache,
    ContactMap& contactMap) const {
  const pandora::ClusterList* pClusterList = NULL;
  PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

  pandora::ClusterList clusterList;
  this->GetRelevantClusterList(*pClusterList, clusterList);

  // Position of each cluster in that list. Contact vectors are built by walking it in order, so an
  // in-place update has to know where in it a cluster belongs.
  ClusterToIndexMap clusterToIndex;
  {
    unsigned int clusterIndex(0);

    for (const pandora::Cluster* const pCluster : clusterList)
      clusterToIndex.emplace(pCluster, clusterIndex++);
  }

  // Create cluster contacts
  pandora::ClusterVector candidateParents, changedClusters;
  pandora::ClusterSet changedClusterSet;

  for (const pandora::Cluster* const pDaughterCluster : clusterList) {
    bool isFullRebuild(true);

    // Identify whether cluster contacts need to be recalculated
    if (!isFirstPass) {
      if (affectedClusters.end() == affectedClusters.find(pDaughterCluster))
        continue;

      // A merge changes exactly two clusters, and a cluster contact is a function of its two clusters
      // and nothing else, so all but a handful of this daughter's contacts are still bit-for-bit what
      // they were. Update those few rather than rebuilding the vector - unless the daughter's own hits
      // moved, in which case every one of its contacts has changed.
      changedClusters = contactCache.GetClustersChangedSince(pDaughterCluster);
      changedClusterSet.clear();
      changedClusterSet.insert(changedClusters.begin(), changedClusters.end());
      isFullRebuild = (changedClusterSet.end() != changedClusterSet.find(pDaughterCluster));

      if (isFullRebuild) {
        const auto pastEntryIter = contactMap.find(pDaughterCluster);

        if (contactMap.end() != pastEntryIter)
          contactMap.erase(pastEntryIter);
      }
    }

    // Apply simple daughter selection cuts
    if ((pDaughterCluster->GetNCaloHits() < m_minDaughterCaloHits) ||
        (pDaughterCluster->GetHadronicEnergy() < m_minDaughterHadronicEnergy) ||
        !this->IsCandidateDaughter(pDaughterCluster)) {
      // A rebuild has already dropped the entry; an update has to, since a rebuild would not put one back.
      if (!isFullRebuild)
        contactMap.erase(pDaughterCluster);

      // No contacts is a state that reflects every merge, so say so and keep the replay list short.
      contactCache.MarkUpToDate(pDaughterCluster);
      continue;
    }

    if (!isFullRebuild) {
      this->UpdateClusterContacts(pDaughterCluster, changedClusters, changedClusterSet, clusterToIndex, contactCache,
                                  contactMap);
      continue;
    }

    // Enumerate the parent candidates that could contribute, in cluster list order, before evaluating any
    // of them. Splitting enumeration from evaluation is what keeps the expensive part off the pairs that
    // cannot pass, and it is the shape a parallel evaluate/serial apply would need.
    const ClusterBoundingBox& daughterBoundingBox(contactCache.GetBoundingBox(pDaughterCluster));
    candidateParents.clear();

    for (const pandora::Cluster* const pParentCluster : clusterList) {
      if (pDaughterCluster == pParentCluster)
        continue;

      if (!this->IsCandidateParent(pDaughterCluster, pParentCluster))
        continue;

      if (!this->CouldPassClusterContactCuts(pDaughterCluster, daughterBoundingBox, pParentCluster, contactCache))
        continue;

      candidateParents.push_back(pParentCluster);
    }

    // Calculate the cluster contact information
    for (const pandora::Cluster* const pParentCluster : candidateParents) {
      const CONTACT contact(this->GetPandora(), pDaughterCluster, pParentCluster, m_contactParameters, contactCache);

      if (this->PassesContactCuts(contact)) {
        contactMap[pDaughterCluster].push_back(contact);
      }
    }

    contactCache.MarkUpToDate(pDaughterCluster);
  }
  isFirstPass = false;

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
void FragmentRemovalBaseAlgorithm<CONTACT>::UpdateClusterContacts(const pandora::Cluster* const pDaughterCluster,
                                                                  const pandora::ClusterVector& changedClusters,
                                                                  const pandora::ClusterSet& changedClusterSet,
                                                                  const ClusterToIndexMap& clusterToIndex,
                                                                  ClusterContactCache& contactCache,
                                                                  ContactMap& contactMap) const {
  const auto mapIter = contactMap.find(pDaughterCluster);

  ContactVector contactVector;

  if (contactMap.end() != mapIter)
    contactVector.swap(mapIter->second);

  // Keep every contact a rebuild would have reproduced unchanged: all of them except those with a cluster
  // some merge has since altered, and those with a parent that has dropped out of the relevant cluster
  // list, which a rebuild would never have reached.
  ContactVector updatedContactVector;
  updatedContactVector.reserve(contactVector.size() + changedClusters.size());

  for (const CONTACT& contact : contactVector) {
    const pandora::Cluster* const pParentCluster(contact.GetParentCluster());

    if (changedClusterSet.end() != changedClusterSet.find(pParentCluster))
      continue;

    if (clusterToIndex.end() == clusterToIndex.find(pParentCluster))
      continue;

    updatedContactVector.push_back(contact);
  }

  // Recompute the contacts those merges invalidated, applying exactly the parent selection and the contact
  // cuts the enumeration loop applies.
  const ClusterBoundingBox& daughterBoundingBox(contactCache.GetBoundingBox(pDaughterCluster));

  for (const pandora::Cluster* const pParentCluster : changedClusters) {
    if ((clusterToIndex.end() == clusterToIndex.find(pParentCluster)) || (pDaughterCluster == pParentCluster) ||
        !this->IsCandidateParent(pDaughterCluster, pParentCluster) ||
        !this->CouldPassClusterContactCuts(pDaughterCluster, daughterBoundingBox, pParentCluster, contactCache)) {
      continue;
    }

    const CONTACT contact(this->GetPandora(), pDaughterCluster, pParentCluster, m_contactParameters, contactCache);

    if (this->PassesContactCuts(contact))
      updatedContactVector.push_back(contact);
  }

  // A rebuild walks the cluster list in order, so the vector has to end up in that order too: the merging
  // candidate search breaks exact ties on the first entry it meets.
  std::sort(updatedContactVector.begin(), updatedContactVector.end(),
            [&clusterToIndex](const CONTACT& lhs, const CONTACT& rhs) {
              return clusterToIndex.at(lhs.GetParentCluster()) < clusterToIndex.at(rhs.GetParentCluster());
            });

  // A rebuild only creates a map entry when it stores a contact, so an emptied vector leaves no entry.
  if (updatedContactVector.empty()) {
    if (contactMap.end() != mapIter)
      contactMap.erase(mapIter);
  } else if (contactMap.end() != mapIter) {
    mapIter->second.swap(updatedContactVector);
  } else {
    contactMap[pDaughterCluster].swap(updatedContactVector);
  }

  contactCache.MarkUpToDate(pDaughterCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
pandora::StatusCode FragmentRemovalBaseAlgorithm<CONTACT>::GetAffectedClusters(
    const ContactMap& contactMap, const pandora::Cluster* const pBestParentCluster,
    const pandora::Cluster* const pBestDaughterCluster, pandora::ClusterSet& affectedClusters) const {
  if (contactMap.end() == contactMap.find(pBestDaughterCluster))
    return pandora::STATUS_CODE_FAILURE;

  affectedClusters.clear();

  for (const auto& mapEntry : contactMap) {
    // Store addresses of all clusters that were in contact with the newly deleted daughter cluster
    if (mapEntry.first == pBestDaughterCluster) {
      for (const CONTACT& contact : mapEntry.second)
        affectedClusters.insert(contact.GetParentCluster());

      continue;
    }

    // Also store addresses of all clusters that contained either the parent or daughter clusters in their own contact
    // vectors
    for (const CONTACT& contact : mapEntry.second) {
      if ((contact.GetParentCluster() == pBestParentCluster) || (contact.GetParentCluster() == pBestDaughterCluster)) {
        affectedClusters.insert(mapEntry.first);
        break;
      }
    }
  }

  return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename CONTACT>
pandora::StatusCode FragmentRemovalBaseAlgorithm<CONTACT>::ReadCommonSettings(const pandora::TiXmlHandle xmlHandle) {
  // Cluster contact parameters
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "ConeCosineHalfAngle1", m_contactParameters.m_coneCosineHalfAngle1));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "CloseHitDistance1", m_contactParameters.m_closeHitDistance1));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "CloseHitDistance2", m_contactParameters.m_closeHitDistance2));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "MinCosOpeningAngle", m_contactParameters.m_minCosOpeningAngle));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "DistanceThreshold", m_contactParameters.m_distanceThreshold));

  // Initial cluster candidate selection
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "MinDaughterCaloHits", m_minDaughterCaloHits));

  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "MinDaughterHadronicEnergy", m_minDaughterHadronicEnergy));

  // Cluster contact cuts
  PANDORA_RETURN_RESULT_IF_AND_IF(
      pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=,
      pandora::XmlHelper::ReadValue(xmlHandle, "ContactCutMaxDistance", m_contactCutMaxDistance));

  return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lc_content

#endif // #ifndef LC_FRAGMENT_REMOVAL_BASE_ALGORITHM_H
