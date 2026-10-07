/**
 *  @file   LCContent/include/LCObjects/LCCaloHit.h
 *
 *  @brief  Header file for the LC calo hit class.
 *
 *  $Log: $
 */
#ifndef LC_CALO_HIT_H
#define LC_CALO_HIT_H 1

#include "Objects/CaloHit.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/ObjectFactory.h"

#include "LCObjects/LCPersistencyConfig.h"

#ifndef LC_PANDORA_FIELD_MAP_PERSISTENCY
#include "Persistency/BinaryFileReader.h"
#include "Persistency/BinaryFileWriter.h"
#include "Persistency/XmlFileReader.h"
#include "Persistency/XmlFileWriter.h"
#endif

namespace lc_content {

namespace calo_hit_bits {

  /// The status bit reserved for the possible beam induced background flag
  constexpr unsigned int POSSIBLE_BIB_BIT = 31u;

  /// The highest status bit available for user defined flags, all bits above this one are reserved for LCContent
  constexpr unsigned int MAX_USER_STATUS_BIT = 30u;
} // namespace calo_hit_bits

namespace utils {

  /**
   *  @brief  Set a single bit of a bitfield to the desired value
   *
   *  @param  bitfield the bitfield to update
   *  @param  bit the bit to set
   *  @param  value the value to set the bit to
   *
   *  @return the updated bitfield
   */
  constexpr unsigned int SetBit(const unsigned int bitfield, const unsigned int bit, const bool value) {
    return (bitfield & ~(1u << bit)) | (static_cast<unsigned int>(value) << bit);
  }

  /**
   *  @brief  Check whether a single bit of a bitfield is set
   *
   *  @param  bitfield the bitfield to check
   *  @param  bit the bit to check
   *
   *  @return whether the bit is set
   */
  constexpr bool CheckBit(const unsigned int bitfield, const unsigned int bit) {
    return ((bitfield >> bit) & 1u) != 0u;
  }
} // namespace utils

/**
 *  @brief  LCCaloHit parameters, the extension point for LC-specific calo hit information
 */
class LCCaloHitParameters : public object_creation::CaloHit::Parameters {
public:
  pandora::InputUInt m_statusBits; ///< The LC-specific calo hit status bits (optional, defaults to 0)
};

/**
 *  @brief  LCCaloHit extension of the CaloHit class for LC-content
 *
 *  Adds a bitfield of status flags on top of the SDK calo hit. The bits above MAX_USER_STATUS_BIT are
 *  reserved for flags defined by LCContent itself, currently only the possible beam induced background
 *  flag. The remaining bits are free for use by downstream clients, which can therefore introduce new
 *  flags without any change to this class.
 */
class LCCaloHit : public object_creation::CaloHit::Object {

public:
  LCCaloHit(const LCCaloHitParameters& parameters)
      : object_creation::CaloHit::Object(parameters),
        m_statusBits(parameters.m_statusBits.IsInitialized() ? parameters.m_statusBits.Get() : 0u) {}

  /**
   *  @brief  Constructor for a calo hit fragment
   *
   *  Forwards to the protected base weighted-copy constructor. Used by LCCaloHitFragmentFactory to
   *  retain the LCCaloHit type when a hit is fragmented.
   *
   *  ATTN: the status bits are inherited from the original calo hit, mirroring the way the SDK calo hit
   *  carries its isolated and possible mip flags over to a fragment.
   *
   *  @param  parameters the calo hit fragment parameters
   */
  LCCaloHit(const object_creation::CaloHitFragment::Parameters& parameters)
      : object_creation::CaloHit::Object(parameters),
        m_statusBits(LCCaloHit::GetOriginalStatusBits(parameters.m_pOriginalCaloHit)) {}

  virtual ~LCCaloHit() = default;

  /**
   *  @brief  Whether the calo hit is possibly due to beam induced background
   *
   *  @return whether the calo hit is a possible beam induced background hit
   */
  bool IsPossibleBIB() const { return utils::CheckBit(m_statusBits, calo_hit_bits::POSSIBLE_BIB_BIT); }

  /**
   *  @brief  Set the possible beam induced background flag
   *
   *  @param  value the value to set the flag to
   */
  void SetPossibleBIB(const bool value = true) {
    m_statusBits = utils::SetBit(m_statusBits, calo_hit_bits::POSSIBLE_BIB_BIT, value);
  }

  /**
   *  @brief  Get a single status bit
   *
   *  @param  bit the bit to check
   *
   *  @return whether the bit is set
   */
  bool GetStatusBit(const unsigned int bit) const {
    if (bit > calo_hit_bits::POSSIBLE_BIB_BIT)
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return utils::CheckBit(m_statusBits, bit);
  }

  /**
   *  @brief  Set a single user defined status bit
   *
   *  ATTN: only the bits up to and including MAX_USER_STATUS_BIT can be set this way, the bits above are
   *  reserved for LCContent and have to be set through their dedicated accessors.
   *
   *  @param  bit the bit to set
   *  @param  value the value to set the bit to
   */
  void SetStatusBit(const unsigned int bit, const bool value = true) {
    if (bit > calo_hit_bits::MAX_USER_STATUS_BIT)
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_statusBits = utils::SetBit(m_statusBits, bit, value);
  }

  /**
   *  @brief  Get the complete status bitfield
   *
   *  @return the status bits
   */
  unsigned int GetStatusBits() const { return m_statusBits; }

private:
  /**
   *  @brief  Get the status bits of a calo hit that may or may not be an LCCaloHit
   *
   *  @param  pCaloHit the address of the calo hit
   *
   *  @return the status bits of the calo hit, or 0 if it is not an LCCaloHit
   */
  static unsigned int GetOriginalStatusBits(const pandora::CaloHit* const pCaloHit) {
    const LCCaloHit* const pLCCaloHit(dynamic_cast<const LCCaloHit*>(pCaloHit));

    return pLCCaloHit ? pLCCaloHit->GetStatusBits() : 0u;
  }

  unsigned int m_statusBits; ///< The LC-specific calo hit status bits
};

/**
 *  @brief  Get the LCCaloHit view of a calo hit
 *
 *  ATTN: returns nullptr rather than throwing, so that content can degrade gracefully when running on
 *  calo hits that were not created through the LCCaloHitFactory.
 *
 *  @param  pCaloHit the address of the calo hit
 *
 *  @return the address of the LCCaloHit, or nullptr if the calo hit is not an LCCaloHit
 */
inline const LCCaloHit* GetLCCaloHit(const pandora::CaloHit* const pCaloHit) {
  return dynamic_cast<const LCCaloHit*>(pCaloHit);
}

/**
 *  @brief  Set the possible beam induced background flag on a calo hit
 *
 *  @param  pCaloHit the address of the calo hit
 *  @param  value the value to set the flag to
 *
 *  @return STATUS_CODE_INVALID_PARAMETER if the calo hit is not an LCCaloHit, STATUS_CODE_SUCCESS otherwise
 */
inline pandora::StatusCode SetPossibleBIB(const pandora::CaloHit* const pCaloHit, const bool value = true) {
  LCCaloHit* const pLCCaloHit(const_cast<LCCaloHit*>(GetLCCaloHit(pCaloHit)));

  if (!pLCCaloHit)
    return pandora::STATUS_CODE_INVALID_PARAMETER;

  pLCCaloHit->SetPossibleBIB(value);

  return pandora::STATUS_CODE_SUCCESS;
}

/**
 *  @brief  LCCaloHitFactory responsible for LCCaloHit creation
 */
class LCCaloHitFactory
    : public pandora::ObjectFactory<object_creation::CaloHit::Parameters, object_creation::CaloHit::Object> {
public:
  /**
   *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
   *
   *  @return the address of the new parameters instance
   */
  Parameters* NewParameters() const { return (new LCCaloHitParameters); }

#ifdef LC_PANDORA_FIELD_MAP_PERSISTENCY
  /**
   *  @brief  Read any additional (derived class only) object parameters from the supplied field map
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  fields the field map, used to extract any additional parameters
   */
  pandora::StatusCode Read(Parameters& parameters, const pandora::FieldMap& fields) const {
    unsigned int statusBits(0u);
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, fields.Get("statusBits", statusBits));

    LCCaloHitParameters& lcCaloHitParameters(dynamic_cast<LCCaloHitParameters&>(parameters));
    lcCaloHitParameters.m_statusBits = statusBits;

    return pandora::STATUS_CODE_SUCCESS;
  }

  /**
   *  @brief  Persist any additional (derived class only) object parameters into the supplied field map
   *
   *  @param  pObject the address of the object to persist
   *  @param  fields the field map to receive the additional parameters
   */
  pandora::StatusCode Write(const Object* const pObject, pandora::FieldMap& fields) const {
    const LCCaloHit* const pLCCaloHit(dynamic_cast<const LCCaloHit*>(pObject));
    if (!pLCCaloHit)
      return pandora::STATUS_CODE_INVALID_PARAMETER;

    fields.Set("statusBits", pLCCaloHit->GetStatusBits());

    return pandora::STATUS_CODE_SUCCESS;
  }
#else
  /**
   *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  fileReader the file reader, used to extract any additional parameters from file
   */
  pandora::StatusCode Read(Parameters& parameters, pandora::FileReader& fileReader) const {
    unsigned int statusBits(0u);

    if (pandora::BINARY == fileReader.GetFileType()) {
      pandora::BinaryFileReader& binaryFileReader(dynamic_cast<pandora::BinaryFileReader&>(fileReader));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(statusBits));
    } else if (pandora::XML == fileReader.GetFileType()) {
      pandora::XmlFileReader& xmlFileReader(dynamic_cast<pandora::XmlFileReader&>(fileReader));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("StatusBits", statusBits));
    } else {
      return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LCCaloHitParameters& lcCaloHitParameters(dynamic_cast<LCCaloHitParameters&>(parameters));
    lcCaloHitParameters.m_statusBits = statusBits;

    return pandora::STATUS_CODE_SUCCESS;
  }

  /**
   *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
   *
   *  @param  pObject the address of the object to persist
   *  @param  fileWriter the file writer
   */
  pandora::StatusCode Write(const Object* const pObject, pandora::FileWriter& fileWriter) const {
    const LCCaloHit* const pLCCaloHit(dynamic_cast<const LCCaloHit*>(pObject));

    if (!pLCCaloHit)
      return pandora::STATUS_CODE_INVALID_PARAMETER;

    const unsigned int statusBits(pLCCaloHit->GetStatusBits());

    if (pandora::BINARY == fileWriter.GetFileType()) {
      pandora::BinaryFileWriter& binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter&>(fileWriter));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(statusBits));
    } else if (pandora::XML == fileWriter.GetFileType()) {
      pandora::XmlFileWriter& xmlFileWriter(dynamic_cast<pandora::XmlFileWriter&>(fileWriter));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("StatusBits", statusBits));
    } else {
      return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    return pandora::STATUS_CODE_SUCCESS;
  }
#endif

  /**
   *  @brief  Create an object with the given parameters
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  pObject to receive the address of the object created
   */
  pandora::StatusCode Create(const Parameters& parameters, const Object*& pObject) const {
    const LCCaloHitParameters& lcCaloHitParameters(dynamic_cast<const LCCaloHitParameters&>(parameters));
    pObject = new LCCaloHit(lcCaloHitParameters);

    return pandora::STATUS_CODE_SUCCESS;
  }
};

/**
 *  @brief  LCCaloHitFragmentFactory responsible for retaining the LCCaloHit type across calo hit fragmentation
 */
class LCCaloHitFragmentFactory : public pandora::ObjectFactory<object_creation::CaloHitFragment::Parameters,
                                                               object_creation::CaloHitFragment::Object> {
public:
  /**
   *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
   *
   *  @return the address of the new parameters instance
   */
  Parameters* NewParameters() const { return (new object_creation::CaloHitFragment::Parameters); }

#ifdef LC_PANDORA_FIELD_MAP_PERSISTENCY
  /**
   *  @brief  Read any additional (derived class only) object parameters from the supplied field map
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  fields the field map, used to extract any additional parameters
   */
  pandora::StatusCode Read(Parameters&, const pandora::FieldMap&) const { return pandora::STATUS_CODE_SUCCESS; }

  /**
   *  @brief  Persist any additional (derived class only) object parameters into the supplied field map
   *
   *  @param  pObject the address of the object to persist
   *  @param  fields the field map to receive the additional parameters
   */
  pandora::StatusCode Write(const Object* const, pandora::FieldMap&) const { return pandora::STATUS_CODE_SUCCESS; }
#else
  /**
   *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  fileReader the file reader, used to extract any additional parameters from file
   */
  pandora::StatusCode Read(Parameters&, pandora::FileReader&) const { return pandora::STATUS_CODE_SUCCESS; }

  /**
   *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
   *
   *  @param  pObject the address of the object to persist
   *  @param  fileWriter the file writer
   */
  pandora::StatusCode Write(const Object* const, pandora::FileWriter&) const { return pandora::STATUS_CODE_SUCCESS; }
#endif

  /**
   *  @brief  Create an object with the given parameters
   *
   *  @param  parameters the parameters to pass in constructor
   *  @param  pObject to receive the address of the object created
   */
  pandora::StatusCode Create(const Parameters& parameters, const Object*& pObject) const {
    pObject = new LCCaloHit(parameters);

    return pandora::STATUS_CODE_SUCCESS;
  }
};

} // namespace lc_content

#endif // #ifndef LC_CALO_HIT_H
