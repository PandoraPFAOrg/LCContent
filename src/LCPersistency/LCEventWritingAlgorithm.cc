/**
 *  @file   LCContent/src/LCPersistency/LCEventWritingAlgorithm.cc
 *
 *  @brief  Implementation of the LCContent event writing algorithm class.
 *
 *  $Log: $
 */

#include "LCPersistency/LCEventWritingAlgorithm.h"

#include "LCObjects/LCTrack.h"

using namespace pandora;

namespace lc_content {

LCEventWritingAlgorithm::LCEventWritingAlgorithm() : EventWritingAlgorithm() {}

//------------------------------------------------------------------------------------------------------------------------------------------

LCEventWritingAlgorithm::~LCEventWritingAlgorithm() {}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LCEventWritingAlgorithm::Initialize() {
  const StatusCode statusCode(this->EventWritingAlgorithm::Initialize());

  if (STATUS_CODE_SUCCESS != statusCode)
    return statusCode;

  m_pEventFileWriter->SetFactory(new LCTrackFactory);

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lc_content
