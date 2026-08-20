/**
 *  @file   LCContent/include/LCPersistency/LCEventReadingAlgorithm.h
 *
 *  @brief  Header file for the LCContent event reading algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CONTENT_EVENT_READING_ALGORITHM_H
#define LC_CONTENT_EVENT_READING_ALGORITHM_H 1

#include "Persistency/EventReadingAlgorithm.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lc_content {

/**
 *  @brief  LCEventReadingAlgorithm class
 */
class LCEventReadingAlgorithm : public EventReadingAlgorithm {
public:
  /**
   *  @brief  Default constructor
   */
  LCEventReadingAlgorithm();

  /**
   *  @brief  Destructor
   */
  ~LCEventReadingAlgorithm();

protected:
  pandora::StatusCode Initialize();
};

} // namespace lc_content

#endif // #ifndef LC_CONTENT_EVENT_READING_ALGORITHM_H
