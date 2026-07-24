/**
 *  @file   LCContent/include/LCTopologicalAssociation/IsolatedHitMergingAlgorithm.h
 * 
 *  @brief  Header file for the isolated hit merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LC_ISOLATED_HIT_MERGING_ALGORITHM_H
#define LC_ISOLATED_HIT_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <vector>

namespace lc_content
{

/**
 *  @brief  IsolatedHitMergingAlgorithm class
 */
class IsolatedHitMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    IsolatedHitMergingAlgorithm();

    /**
     *  @brief  Per-cluster cached quantities to avoid recomputing inside the hit loop.
     */
    struct ClusterCache
    {
        const pandora::Cluster         *pCluster;    ///< Pointer to the cluster (nullptr if deleted)
        pandora::CartesianVector        centroid;    ///< Centroid at inner pseudo-layer
        float                           energy;      ///< Cluster energy (raw or plugin-corrected hadronic)
        pandora::CartesianVector        direction;   ///< Cluster initial direction (for opening-angle preselection)
        bool                            isEm;        ///< Whether the cluster is identified as an EM shower
        bool                            isCharged;   ///< Whether the cluster has an associated track

        ClusterCache()
            : pCluster(nullptr)
            , centroid(0.f, 0.f, 0.f)
            , energy(0.f)
            , direction(0.f, 0.f, 1.f)
            , isEm(false)
            , isCharged(false)
        {}
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Build a cache of per-cluster quantities, optionally filtering photons.
     *
     *  @param  clusterVector   input cluster vector
     *  @param  clusterCaches   output vector of ClusterCache entries
     */
    pandora::StatusCode BuildClusterCache(const pandora::ClusterVector &clusterVector,
                                          std::vector<ClusterCache> &clusterCaches) const;

    /**
     *  @brief  Get distance between a calo hit and a cluster using precomputed cache.
     *
     *  @param  cache        precomputed cluster info
     *  @param  pCaloHit     address of the calo hit
     *  @param  distance     output: distance metric (angular 1-cosA or Cartesian mm)
     */
    pandora::StatusCode GetDistanceToHit(const ClusterCache &cache,
                                         const pandora::CaloHit *const pCaloHit,
                                         float &distance) const;

    /**
     *  @brief  Get the energy of a cluster: raw hadronic energy, or the EnergyCorrections-plugin
     *          corrected hadronic energy when UseCorrectedHadronicEnergy is set.
     *
     *  @param  pCluster  address of the cluster
     *  @param  energy    output energy
     */
    pandora::StatusCode GetClusterEnergy(const pandora::Cluster *const pCluster, float &energy) const;

    bool                    m_shouldUseCurrentClusterList;  ///< Whether to use clusters from the current list in the algorithm
    pandora::StringVector   m_additionalClusterListNames;   ///< Additional cluster lists from which to consider clusters

    unsigned int            m_minHitsInCluster;             ///< Min number of hits allowed in a cluster - smaller clusters will be split up
    float                   m_maxRecombinationDistance;     ///< Max distance between calo hit and cluster to allow addition of hit
    float                   m_minCosOpeningAngle;           ///< Min cos(angle) between hit and cluster directions to allow addition of hit

    bool                    m_useCorrectedHadronicEnergy;   ///< Use the plugin-corrected hadronic energy instead of raw GetHadronicEnergy() (default: false)
    bool                    m_ignorePhotons;                ///< Whether to ignore photons when merging isolated hits (default: false)
    bool                    m_ignoreCharged;                ///< Whether to ignore charged clusters when merging isolated hits (default: false)
    bool                    m_useAngularDistance;           ///< Whether to use angular (1-cos(angle)) distance instead of Cartesian (default: false)
};

} // namespace lc_content

#endif // #ifndef LC_ISOLATED_HIT_MERGING_ALGORITHM_H
