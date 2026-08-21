/**
 *  @file   LCContent/include/LCPersistency/LCEventWritingAlgorithm.h
 *
 *  @brief  Header file for the LCContent event writing algorithm class.
 *
 *  $Log: $
 */
#ifndef LC_CONTENT_EVENT_WRITING_ALGORITHM_H
#define LC_CONTENT_EVENT_WRITING_ALGORITHM_H 1

#include "Persistency/EventWritingAlgorithm.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lc_content {

/**
 *  @brief  LCEventWritingAlgorithm class
 */
class LCEventWritingAlgorithm : public EventWritingAlgorithm {
public:
  /**
   *  @brief  Default constructor
   */
  LCEventWritingAlgorithm();

  /**
   *  @brief  Destructor
   */
  ~LCEventWritingAlgorithm();

protected:
  pandora::StatusCode Initialize();
};

} // namespace lc_content

#endif // #ifndef LC_CONTENT_EVENT_WRITING_ALGORITHM_H
