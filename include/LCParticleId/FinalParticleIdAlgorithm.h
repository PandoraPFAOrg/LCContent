/**
 *  @file   LCContent/include/LCParticleId/FinalParticleIdAlgorithm.h
 *
 *  @brief  Header file for the final particle id algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_FINAL_PARTICLE_ID_ALGORITHM_H
#define LC_FINAL_PARTICLE_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lc_content {

/**
 *  @brief  FinalParticleIdAlgorithm class
 */
class FinalParticleIdAlgorithm : public pandora::Algorithm {
private:
  pandora::StatusCode Run();
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  /**
   *  @brief  Determine whether any track in a specified tracklist has a sibling
   *
   *  @param  trackList the tracklist
   *
   *  @return boolean
   */
  bool ContainsSiblingTrack(const pandora::TrackList& trackList) const;
};

} // namespace lc_content

#endif // #ifndef LC_FINAL_PARTICLE_ID_ALGORITHM_H
