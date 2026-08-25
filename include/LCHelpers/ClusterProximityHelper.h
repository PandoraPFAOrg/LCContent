/**
 *  @file   LCContent/include/LCHelpers/ClusterProximityHelper.h
 *
 *  @brief  Header file for the cluster proximity helper classes.
 *
 *  $Log: $
 */
#ifndef LC_CLUSTER_PROXIMITY_HELPER_H
#define LC_CLUSTER_PROXIMITY_HELPER_H 1

#include "Objects/CartesianVector.h"

#include "Pandora/PandoraInternal.h"

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

namespace lc_content {

/**
 *  @brief  ClusterBoundingBox class, the axis-aligned box enclosing the calo hit positions of a cluster
 *
 *          Several algorithms compare every pair of clusters with a doubly-nested loop over the hits of
 *          each, then throw the answer away because the two clusters are nowhere near one another. A
 *          bounding box answers "nowhere near" in constant time.
 *
 *          Every separation the box reports is a lower bound on the corresponding hit-hit distance: a
 *          box reported as separated really is separated, while a box not reported as separated may
 *          still be. A caller must read false as "not certain" and fall through to the hit loop.
 */
class ClusterBoundingBox {
public:
  /**
   *  @brief  Default constructor, giving an uninitialized box that is never separated from anything
   */
  ClusterBoundingBox();

  /**
   *  @brief  Constructor, enclosing the calo hits of the provided cluster
   *
   *  @param  pCluster address of the cluster
   */
  explicit ClusterBoundingBox(const pandora::Cluster* const pCluster);

  /**
   *  @brief  Constructor, enclosing the provided calo hits
   *
   *  @param  caloHitList the calo hits
   */
  explicit ClusterBoundingBox(const pandora::CaloHitList& caloHitList);

  /**
   *  @brief  Whether every hit in this box is further than the provided distance from every hit in another box
   *
   *  @param  rhs the second box
   *  @param  distance the distance
   *
   *  @return boolean, false whenever the answer is not certain
   */
  bool IsSeparatedFrom(const ClusterBoundingBox& rhs, const float distance) const;

  /**
   *  @brief  Whether every hit in this box is further than the square root of the provided distance from a point
   *
   *          The squared form is what a hit loop already has to hand, so it saves the caller a square root
   *          in its innermost test.
   *
   *  @param  point the point
   *  @param  distanceSquared the squared distance
   *
   *  @return boolean, false whenever the answer is not certain
   */
  bool IsSeparatedFromSquared(const pandora::CartesianVector& point, const float distanceSquared) const;

private:
  /**
   *  @brief  Get the squared distance from a point to the box, zero if the point lies inside it
   *
   *  @param  point the point
   *
   *  @return the squared separation
   */
  float GetSeparationSquared(const pandora::CartesianVector& point) const;

  /**
   *  @brief  Get the squared distance between two boxes, zero if they overlap
   *
   *  @param  rhs the second box
   *
   *  @return the squared separation
   */
  float GetSeparationSquared(const ClusterBoundingBox& rhs) const;

  /**
   *  @brief  Grow the box to enclose the provided calo hits
   *
   *  @param  caloHitList the calo hits
   */
  void Enclose(const pandora::CaloHitList& caloHitList);

  /**
   *  @brief  Whether a squared lower bound is certainly beyond a squared threshold
   *
   *  @param  separationSquared the squared lower bound
   *  @param  distanceSquared the squared threshold
   *
   *  @return boolean
   */
  static bool IsCertainlyBeyond(const float separationSquared, const float distanceSquared);

  pandora::CartesianVector m_min; ///< Lower corner of the box
  pandora::CartesianVector m_max; ///< Upper corner of the box
  bool m_isInitialized;           ///< Whether the box encloses at least one hit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Position of each cluster in the list an algorithm is iterating
 *
 *          Contact vectors are built by walking that list in order, so anything that updates one in place
 *          rather than rebuilding it needs to be able to ask where a cluster sits in it.
 */
using ClusterToIndexMap = std::unordered_map<const pandora::Cluster*, unsigned int>;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  One bounding box per occupied pseudo layer of a cluster, in pseudo layer order
 */
using ClusterLayerBoundingBoxVector = std::vector<ClusterBoundingBox>;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ClusterContactCache class, the state the fragment removal algorithms carry from one pass to the next
 *
 *          Those algorithms merge one cluster pair per pass and then reconsider the clusters the merge
 *          disturbed. Two things make that cheap, and both need memory of what has happened so far.
 *
 *          Bounding boxes are the first. Nothing inside a fragment removal Run creates a cluster, so the
 *          list only ever shrinks and cluster addresses are never recycled; a box stays valid until its
 *          cluster gains hits, which happens only to the parent of a merge.
 *
 *          The merge log is the second, and it is what makes an incremental contact update safe. A cluster
 *          contact is a function of its two clusters and nothing else, so a merge invalidates only the
 *          contacts involving the two clusters it touched. A daughter that is brought up to date every pass
 *          therefore needs one contact recomputed rather than a whole vector. The catch is that a daughter
 *          the algorithm skipped for several passes has missed several merges, and a full rebuild would
 *          have quietly repaired all of them at once. Recording every merge, and giving each cluster a
 *          watermark saying how many of them its contacts already reflect, lets the update replay exactly
 *          the ones it missed - which is what makes it produce the same answer a rebuild would, and never
 *          more work than a rebuild would have cost.
 *
 *          Callers hold on to the references the two box getters return while asking for other clusters'
 *          boxes, which is safe because unordered_map keeps references to its elements valid across a
 *          rehash. A container without that guarantee cannot be substituted here.
 */
class ClusterContactCache {
public:
  /**
   *  @brief  Get the bounding box of a cluster, computing it on first use
   *
   *  @param  pCluster address of the cluster
   *
   *  @return the bounding box
   */
  const ClusterBoundingBox& GetBoundingBox(const pandora::Cluster* const pCluster);

  /**
   *  @brief  Get one bounding box per occupied pseudo layer of a cluster, in pseudo layer order, computing them on
   * first use
   *
   *          A pseudo layer is a thin shell, so these bound a cluster's hits far more tightly than one box
   *          around the whole cluster does. They depend on the cluster alone, and a cluster is the parent of
   *          many pairs, so they are worth keeping rather than rebuilding per pair.
   *
   *          The vector is built by walking the cluster's ordered calo hit list, so its entries pair up with
   *          that list in order. A caller stepping through the two together is what makes the pairing correct.
   *
   *  @param  pCluster address of the cluster
   *
   *  @return the per pseudo layer bounding boxes
   */
  const ClusterLayerBoundingBoxVector& GetLayerBoundingBoxes(const pandora::Cluster* const pCluster);

  /**
   *  @brief  Record a merge: the parent has gained the daughter's hits and the daughter has been deleted
   *
   *  @param  pParentCluster address of the parent cluster
   *  @param  pDaughterCluster address of the daughter cluster, no longer dereferenceable after this call
   */
  void RecordMerge(const pandora::Cluster* const pParentCluster, const pandora::Cluster* const pDaughterCluster);

  /**
   *  @brief  Get the number of merges recorded so far
   *
   *  @return the number of merges
   */
  unsigned int GetNMerges() const;

  /**
   *  @brief  Get every cluster changed by the merges a cluster's contacts have not yet seen
   *
   *          Returned in merge order and without repeats. These are exactly the parents whose contacts
   *          with the provided cluster have to be recomputed; every other contact it holds is unchanged.
   *
   *  @param  pCluster address of the cluster
   *
   *  @return the changed clusters
   */
  pandora::ClusterVector GetClustersChangedSince(const pandora::Cluster* const pCluster) const;

  /**
   *  @brief  Declare a cluster's contacts to reflect every merge recorded so far
   *
   *  @param  pCluster address of the cluster
   */
  void MarkUpToDate(const pandora::Cluster* const pCluster);

private:
  using ClusterToBoundingBoxMap = std::unordered_map<const pandora::Cluster*, ClusterBoundingBox>;
  using ClusterToLayerBoundingBoxMap = std::unordered_map<const pandora::Cluster*, ClusterLayerBoundingBoxVector>;
  using ClusterToWatermarkMap = std::unordered_map<const pandora::Cluster*, unsigned int>;
  using ClusterMerge = std::pair<const pandora::Cluster*, const pandora::Cluster*>;

  ClusterToBoundingBoxMap m_boundingBoxes;           ///< The cached whole-cluster bounding boxes
  ClusterToLayerBoundingBoxMap m_layerBoundingBoxes; ///< The cached per pseudo layer bounding boxes
  ClusterToWatermarkMap m_watermarks;                ///< Number of merges each cluster's cached contacts reflect
  std::vector<ClusterMerge> m_merges;                ///< The parent and daughter of every merge, in order
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterBoundingBox::IsCertainlyBeyond(const float separationSquared, const float distanceSquared) {
  // The margin makes the lower-bound guarantee survive floating point rounding: the box separation and the
  // hit-hit distance it bounds are formed by different expressions, and a compiler is free to contract
  // a*a + b into an fma in one and not the other, which can move a result by an ulp. The margin is many ulp
  // wide and, at the distances these cuts use, a few microns of slack, so a pair sitting on the threshold is
  // always evaluated rather than skipped.
  const float boundSafetyFactor(1.f - 1e-5f);
  return ((separationSquared * boundSafetyFactor) > distanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterBoundingBox::IsSeparatedFromSquared(const pandora::CartesianVector& point,
                                                       const float distanceSquared) const {
  if (!m_isInitialized)
    return false;

  return ClusterBoundingBox::IsCertainlyBeyond(this->GetSeparationSquared(point), distanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ClusterBoundingBox::IsSeparatedFrom(const ClusterBoundingBox& rhs, const float distance) const {
  if (!m_isInitialized || !rhs.m_isInitialized)
    return false;

  return ClusterBoundingBox::IsCertainlyBeyond(this->GetSeparationSquared(rhs), distance * distance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ClusterBoundingBox::GetSeparationSquared(const pandora::CartesianVector& point) const {
  const float dx(std::max(0.f, std::max(m_min.GetX() - point.GetX(), point.GetX() - m_max.GetX())));
  const float dy(std::max(0.f, std::max(m_min.GetY() - point.GetY(), point.GetY() - m_max.GetY())));
  const float dz(std::max(0.f, std::max(m_min.GetZ() - point.GetZ(), point.GetZ() - m_max.GetZ())));

  return ((dx * dx) + (dy * dy) + (dz * dz));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ClusterBoundingBox::GetSeparationSquared(const ClusterBoundingBox& rhs) const {
  const float dx(std::max(0.f, std::max(m_min.GetX() - rhs.m_max.GetX(), rhs.m_min.GetX() - m_max.GetX())));
  const float dy(std::max(0.f, std::max(m_min.GetY() - rhs.m_max.GetY(), rhs.m_min.GetY() - m_max.GetY())));
  const float dz(std::max(0.f, std::max(m_min.GetZ() - rhs.m_max.GetZ(), rhs.m_min.GetZ() - m_max.GetZ())));

  return ((dx * dx) + (dy * dy) + (dz * dz));
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ClusterContactCache::GetNMerges() const { return static_cast<unsigned int>(m_merges.size()); }

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ClusterContactCache::MarkUpToDate(const pandora::Cluster* const pCluster) {
  m_watermarks[pCluster] = this->GetNMerges();
}

} // namespace lc_content

#endif // #ifndef LC_CLUSTER_PROXIMITY_HELPER_H
