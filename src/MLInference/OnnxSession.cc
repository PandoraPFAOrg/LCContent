/**
 *  @file   LCContent/src/MLInference/OnnxSession.cc
 *
 *  @brief  Implementation of OnnxSession -- the ONLY translation unit that includes onnxruntime.
 */
#include "MLInference/OnnxSession.h"

#include <onnxruntime_cxx_api.h>

#include <iostream>

namespace lc_content {

struct OnnxSession::Impl {
  Ort::Env env{ORT_LOGGING_LEVEL_WARNING, "OnnxSession"};
  Ort::MemoryInfo mem{Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault)};
  std::unique_ptr<Ort::Session> session;
  std::vector<std::string> inputNames, outputNames;
  bool valid = false;
};

//------------------------------------------------------------------------------------------------------------------------------------------

OnnxSession::OnnxSession(const std::string &modelPath, const int intraOpThreads) : m_impl(new Impl()) {
  try {
    Ort::SessionOptions options;
    options.SetIntraOpNumThreads(intraOpThreads);
    m_impl->session.reset(new Ort::Session(m_impl->env, modelPath.c_str(), options)); // throws on bad path

    Ort::AllocatorWithDefaultOptions allocator;
    for (std::size_t i = 0; i < m_impl->session->GetInputCount(); ++i)
      m_impl->inputNames.emplace_back(m_impl->session->GetInputNameAllocated(i, allocator).get());
    for (std::size_t i = 0; i < m_impl->session->GetOutputCount(); ++i)
      m_impl->outputNames.emplace_back(m_impl->session->GetOutputNameAllocated(i, allocator).get());

    m_impl->valid = true;
  } catch (const Ort::Exception &exception) {
    std::cout << "OnnxSession: failed to load '" << modelPath << "' (" << exception.what() << ")." << std::endl;
    m_impl->valid = false;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

OnnxSession::~OnnxSession() = default;

//------------------------------------------------------------------------------------------------------------------------------------------

bool OnnxSession::IsValid() const { return m_impl->valid; }

std::size_t OnnxSession::InputCount() const { return m_impl->inputNames.size(); }

std::size_t OnnxSession::OutputCount() const { return m_impl->outputNames.size(); }

const std::vector<std::string> &OnnxSession::InputNames() const { return m_impl->inputNames; }

const std::vector<std::string> &OnnxSession::OutputNames() const { return m_impl->outputNames; }

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::int64_t> OnnxSession::InputShape(const std::size_t index) const {
  if (!m_impl->valid || index >= m_impl->session->GetInputCount())
    return {};

  return m_impl->session->GetInputTypeInfo(index).GetTensorTypeAndShapeInfo().GetShape();
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::vector<float>> OnnxSession::Run(const std::vector<Input> &inputs) const {
  if (!m_impl->valid)
    return {};

  // ONNX Runtime's CreateTensor takes a non-const data pointer but does not modify inputs; the
  // const_cast is safe (the caller owns const data for the duration of the call).
  std::vector<Ort::Value> ortInputs;
  ortInputs.reserve(inputs.size());
  for (const Input &input : inputs) {
    std::size_t nElements = 1;
    for (const std::int64_t dim : input.shape)
      nElements *= static_cast<std::size_t>(dim);

    switch (input.type) {
    case Input::Type::Float:
      ortInputs.emplace_back(Ort::Value::CreateTensor<float>(m_impl->mem, const_cast<float *>(static_cast<const float *>(input.data)),
                                                             nElements, input.shape.data(), input.shape.size()));
      break;
    case Input::Type::Int64:
      ortInputs.emplace_back(
          Ort::Value::CreateTensor<std::int64_t>(m_impl->mem, const_cast<std::int64_t *>(static_cast<const std::int64_t *>(input.data)),
                                                 nElements, input.shape.data(), input.shape.size()));
      break;
    case Input::Type::Bool:
      ortInputs.emplace_back(Ort::Value::CreateTensor<bool>(m_impl->mem, const_cast<bool *>(static_cast<const bool *>(input.data)),
                                                            nElements, input.shape.data(), input.shape.size()));
      break;
    }
  }

  std::vector<const char *> inputNamesC, outputNamesC;
  for (const std::string &name : m_impl->inputNames)
    inputNamesC.push_back(name.c_str());
  for (const std::string &name : m_impl->outputNames)
    outputNamesC.push_back(name.c_str());

  try {
    auto outputs = m_impl->session->Run(Ort::RunOptions{nullptr}, inputNamesC.data(), ortInputs.data(), ortInputs.size(),
                                        outputNamesC.data(), outputNamesC.size());

    std::vector<std::vector<float>> result;
    result.reserve(outputs.size());
    for (Ort::Value &output : outputs) {
      const std::size_t count = output.GetTensorTypeAndShapeInfo().GetElementCount();
      const float *const data = output.GetTensorMutableData<float>();
      result.emplace_back(data, data + count);
    }
    return result;
  } catch (const Ort::Exception &exception) {
    std::cout << "OnnxSession: inference failed (" << exception.what() << ")." << std::endl;
    return {};
  }
}

} // namespace lc_content
