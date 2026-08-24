/**
 *  @file   LCContent/include/MLInference/OnnxSession.h
 *
 *  @brief  Thin PIMPL wrapper around an ONNX Runtime session.  The onnxruntime dependency is
 *          confined to OnnxSession.cc, so the ML-inference algorithms (and their headers) never
 *          include any onnxruntime header.
 */
#ifndef LC_MLINFERENCE_ONNX_SESSION_H
#define LC_MLINFERENCE_ONNX_SESSION_H 1

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace lc_content {

/**
 *  @brief  OnnxSession class: owns an ONNX Runtime session and runs inference on flat buffers.
 */
class OnnxSession {
public:
  /**
   *  @brief  Load the model at modelPath.  On failure the session is left INVALID (see IsValid());
   *          the constructor never throws.
   *
   *  @param  modelPath        path to the .onnx model file
   *  @param  intraOpThreads   ONNX Runtime intra-op thread count (default 1)
   */
  explicit OnnxSession(const std::string& modelPath, int intraOpThreads = 1);

  ~OnnxSession();

  OnnxSession(const OnnxSession&) = delete;
  OnnxSession& operator=(const OnnxSession&) = delete;

  /**
   *  @brief  Whether the model loaded successfully.
   */
  bool IsValid() const;

  std::size_t InputCount() const;
  std::size_t OutputCount() const;
  const std::vector<std::string>& InputNames() const;
  const std::vector<std::string>& OutputNames() const;

  /**
   *  @brief  The model's declared shape for input `index` (e.g. to validate a fixed token dimension).
   *          Empty if the session is invalid or the index is out of range.
   */
  std::vector<std::int64_t> InputShape(std::size_t index) const;

  /**
   *  @brief  One model input: its row-major shape, element type, and a pointer to caller-owned data
   *          which must stay valid for the duration of the Run() call.
   */
  struct Input {
    enum class Type { Float, Int64, Bool };

    std::vector<std::int64_t> shape;
    Type type = Type::Float;
    const void* data = nullptr;
  };

  /**
   *  @brief  Run inference.  Inputs are positional (inputs[i] feeds InputNames()[i]).  Each model
   *          output is returned flattened to float (these models' outputs are float).
   *
   *  @param  inputs  the model inputs, in the model's input order
   *
   *  @return one float vector per model output; empty on failure (which is logged)
   */
  std::vector<std::vector<float>> Run(const std::vector<Input>& inputs) const;

private:
  struct Impl; ///< holds the Ort::Env / Session / MemoryInfo / I/O names
  std::unique_ptr<Impl> m_impl;
};

} // namespace lc_content

#endif // #ifndef LC_MLINFERENCE_ONNX_SESSION_H
