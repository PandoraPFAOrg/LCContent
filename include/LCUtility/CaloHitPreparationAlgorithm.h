/**
 *  @file   LCContent/include/LCUtility/CaloHitPreparationAlgorithm.h
 *
 *  @brief  Header file for the calo hit preparation algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CALO_HIT_PREPARATION_ALGORITHM_H
#define LC_CALO_HIT_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <string>

namespace lc_content {

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CaloHitPreparationAlgorithm class
 */
class CaloHitPreparationAlgorithm : public pandora::Algorithm {
public:
  typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 4> HitKDTree4D;
  typedef KDTreeNodeInfoT<const pandora::CaloHit*, 4> HitKDNode4D;

  /**
   *  @brief Default constructor
   */
  CaloHitPreparationAlgorithm();

  /**
   * @brief Destructor
   */
  ~CaloHitPreparationAlgorithm();

private:
  pandora::StatusCode Run();

  /**
   *  @brief  Initialize a kd-tree of the input hits to the preparation alg.
   *
   *  @param  pCaloHitList -- the calorimeter hit list
   */
  void InitializeKDTree(const pandora::CaloHitList* const pCaloHitList);

  /**
   *  @brief  Calculate calo hit properties for a particular calo hit, through comparison with an ordered list of other
   * hits.
   *
   *  @param  pCaloHit the calo hit
   *  @param  pOrderedCaloHitList the ordered calo hit list
   */
  void CalculateCaloHitProperties(const pandora::CaloHit* const pCaloHit,
                                  const pandora::OrderedCaloHitList& orderedCaloHitList);

  /**
   *  @brief  Count number of "nearby" hits using the isolation scheme
   *
   *  @param  searchLayer the pseudolayer to search in
   *  @param  pCaloHit the calo hit
   *
   *  @return the number of nearby hits
   */
  unsigned int IsolationCountNearbyHits(unsigned int searchLayer, const pandora::CaloHit* const pCaloHit);

  /**
   *  @brief  Count number of "nearby" hits using the mip identification scheme
   *
   *  @param  searchLayer the pseudolayer to search in
   *  @param  pCaloHit the calo hit
   *
   *  @return the number of nearby hits
   */
  unsigned int MipCountNearbyHits(unsigned int searchLayer, const pandora::CaloHit* const pCaloHit);

  /**
   *  @brief  Sum electromagnetic energy of "nearby" hits (excluding the hit itself), using the
   *          same transverse-distance scheme as the isolation calculation
   *
   *  @param  searchLayer the pseudolayer to search in
   *  @param  pCaloHit the calo hit
   *
   *  @return the summed nearby electromagnetic energy, units GeV
   */
  float BIBSumNearbyEnergy(unsigned int searchLayer, const pandora::CaloHit* const pCaloHit);

  /**
   *  @brief  Get the BIB energy density cut for a calo hit, parametrized vs pseudolayer and
   *          selected by subdetector (ECAL/HCAL) and region (barrel/endcap). A subdetector/region
   *          with no cuts file configured, a hit type with no cut table (e.g. muon hits), or a
   *          pseudolayer beyond the end of the table all give 0, i.e. the hit is never flagged.
   *
   *  @param  pCaloHit the calo hit
   *
   *  @return the energy density cut, units GeV/mm^2/layer
   */
  float GetBIBEnergyDensityCut(const pandora::CaloHit* const pCaloHit) const;

  /**
   *  @brief  Load per-pseudolayer BIB energy density cuts from a text file with lines of
   *          "pseudoLayer cutEnergyDensity" ('#' starts a comment); layers not listed get 0.
   *          An empty file name leaves the cuts empty (that subdetector/region is disabled).
   *
   *  @param  fileName the cuts file name
   *  @param  cuts the vector to fill, indexed by pseudolayer
   */
  pandora::StatusCode LoadBIBCuts(const std::string& fileName, pandora::FloatVector& cuts) const;

  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  float m_caloHitMaxSeparation2; ///< Max separation to consider associations between hits, units mm (used squared)
  float m_isolationCaloHitMaxSeparation2; ///< Max separation considered when identifying isolated hits, units mm (used
                                          ///< squared)

  unsigned int m_isolationNLayers;       ///< Number of adjacent layers to use in isolation calculation
  float m_isolationCutDistanceFine2;     ///< Fine granularity isolation cut distance, units mm (used squared)
  float m_isolationCutDistanceCoarse2;   ///< Coarse granularity isolation cut distance, units mm (used squared)
  float m_isolationSearchSafetyFactor;   ///< Safety factor, applied to isolation cut distance, to define kd-tree search
                                         ///< region
  unsigned int m_isolationMaxNearbyHits; ///< Max number of "nearby" hits for a hit to be considered isolated

  float m_mipLikeMipCut;                ///< Mip equivalent energy cut for hit to be flagged as possible mip
  unsigned int m_mipNCellsForNearbyHit; ///< Separation (in calo cells) for hits to be declared "nearby"
  unsigned int m_mipMaxNearbyHits;      ///< Max number of "nearby" hits for hit to be flagged as possible mip

  unsigned int m_bibNLayers;              ///< Number of adjacent layers to use in BIB energy density calculation
  float m_bibCaloHitMaxSeparationSquared; ///< Max separation considered when summing nearby BIB energy, stored squared,
                                          ///< units mm^2 (XML value BIBCaloHitMaxSeparation is in mm)
  float m_bibCutDistanceSquared; ///< Transverse distance defining the BIB neighbourhood, stored squared, units mm^2
                                 ///< (XML value BIBCutDistance is in mm)
  float m_bibSearchSafetyFactor; ///< Safety factor, applied to BIB cut distance, to define kd-tree search region

  std::string m_bibCutsFileECalBarrel; ///< Per-pseudolayer BIB cuts file, ECAL barrel (empty = disabled)
  std::string m_bibCutsFileECalEndcap; ///< Per-pseudolayer BIB cuts file, ECAL endcap (empty = disabled)
  std::string m_bibCutsFileHCalBarrel; ///< Per-pseudolayer BIB cuts file, HCAL barrel (empty = disabled)
  std::string m_bibCutsFileHCalEndcap; ///< Per-pseudolayer BIB cuts file, HCAL endcap (empty = disabled)
  pandora::FloatVector
      m_bibCutsECalBarrel; ///< BIB energy density cuts vs pseudolayer, ECAL barrel, units GeV/mm^2/layer
  pandora::FloatVector
      m_bibCutsECalEndcap; ///< BIB energy density cuts vs pseudolayer, ECAL endcap, units GeV/mm^2/layer
  pandora::FloatVector
      m_bibCutsHCalBarrel; ///< BIB energy density cuts vs pseudolayer, HCAL barrel, units GeV/mm^2/layer
  pandora::FloatVector
      m_bibCutsHCalEndcap; ///< BIB energy density cuts vs pseudolayer, HCAL endcap, units GeV/mm^2/layer

  bool m_bibTimingCutEnabled;  ///< Whether to also flag hits with a time outside the BIB timing window
  float m_bibTimingWindowLow;  ///< Lower edge of the accepted hit time window, units ns
  float m_bibTimingWindowHigh; ///< Upper edge of the accepted hit time window, units ns

  std::vector<HitKDNode4D>* m_hitNodes4D; ///< nodes for the KD tree (used for filling)
  HitKDTree4D* m_hitsKdTree4D;            ///< the kd-tree itself, 4D in x,y,z,pseudolayer
};

} // namespace lc_content

#endif // #ifndef LC_CALO_HIT_PREPARATION_ALGORITHM_H
