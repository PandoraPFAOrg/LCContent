/**
 *  @file   LCContent/src/LCPersistency/LCEventReadingAlgorithm.cc
 *
 *  @brief  Implementation of the LCContent event reading algorithm class.
 *
 *  $Log: $
 */

#include "LCPersistency/LCEventReadingAlgorithm.h"

#include "LCObjects/LCTrack.h"

using namespace pandora;

namespace lc_content {

LCEventReadingAlgorithm::LCEventReadingAlgorithm() : EventReadingAlgorithm() {}

//------------------------------------------------------------------------------------------------------------------------------------------

LCEventReadingAlgorithm::~LCEventReadingAlgorithm() {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEventReadingAlgorithm::Initialize() {
  const StatusCode statusCode(this->EventReadingAlgorithm::Initialize());

  if (STATUS_CODE_SUCCESS != statusCode)
    return statusCode;

  m_pEventFileReader->SetFactory(new LCTrackFactory);

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lc_content
