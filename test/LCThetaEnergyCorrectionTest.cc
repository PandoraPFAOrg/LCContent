/**
 *  @file   LCContent/test/LCThetaEnergyCorrectionTest.cc
 *
 *  @brief  Unit tests for the theta-energy binned non-linearity correction lookup.
 *
 *          The lookup registry is process-global and cannot be cleared, so each case
 *          runs in its own process; argv[1] selects the case and test/CMakeLists.txt
 *          registers one CTest entry per case.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Objects/CartesianVector.h"
#include "Pandora/Pandora.h"
#include "Pandora/StatusCodes.h"

#include "LCContent.h"
#include "LCPlugins/LCEnergyCorrectionPlugins.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace pandora;
using lc_content::LCEnergyCorrectionPlugins;

namespace {

int g_failures(0);

void Check(const bool condition, const std::string& what) {
  std::cout << (condition ? "  pass  " : "  FAIL  ") << what << std::endl;
  if (!condition)
    ++g_failures;
}

void CheckClose(const float actual, const float expected, const std::string& what) {
  const bool ok(std::fabs(actual - expected) < 1e-4f * std::max(1.f, std::fabs(expected)));
  std::cout << (ok ? "  pass  " : "  FAIL  ") << what << "  (expected " << expected << ", got " << actual << ")"
            << std::endl;
  if (!ok)
    ++g_failures;
}

CartesianVector DirAtTheta(const float theta) { return CartesianVector(std::sin(theta), 0.f, std::cos(theta)); }

float Factor(const float theta, const float energy) {
  return LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, DirAtTheta(theta), energy) / energy;
}

// 2 theta bins x 3 energy bins, deliberately non-square so a transposed index is caught.
//   theta bins: [0,1) [1,2)        energy bins: [0,10) [10,50) [50,100)
void RegisterAsymmetricTable(const Pandora& pandora) {
  const FloatVector thetaEdges{0.f, 1.f, 2.f};
  const FloatVector energyEdges{0.f, 10.f, 50.f, 100.f};
  const FloatVector factors{1.1f, 1.2f, 1.3f,  // thetaBin 0
                            2.1f, 2.2f, 2.3f}; // thetaBin 1
  if (STATUS_CODE_SUCCESS !=
      LCContent::RegisterNonLinearityEnergyCorrection(pandora, "Asym", HADRONIC, thetaEdges, energyEdges, factors))
    throw StatusCodeException(STATUS_CODE_FAILURE);
}

bool ExpectRejected(const Pandora& pandora, const std::string& name, const FloatVector& thetaEdges,
                    const FloatVector& energyEdges, const FloatVector& factors, const std::string& what) {
  try {
    const StatusCode sc(
        LCContent::RegisterNonLinearityEnergyCorrection(pandora, name, HADRONIC, thetaEdges, energyEdges, factors));
    if (STATUS_CODE_SUCCESS == sc) {
      std::cout << "  FAIL  " << what << "  (accepted a malformed table)" << std::endl;
      ++g_failures;
      return false;
    }
    std::cout << "  pass  " << what << "  (returned " << StatusCodeToString(sc) << ")" << std::endl;
  } catch (const StatusCodeException& e) {
    std::cout << "  pass  " << what << "  (threw " << e.ToString() << ")" << std::endl;
  }
  return true;
}

//----------------------------------------------------------------------------------------------------------------

// A: bin lookup and row-major indexing.
void CaseLookup(const Pandora& pandora) {
  RegisterAsymmetricTable(pandora);

  CheckClose(Factor(0.5f, 5.f), 1.1f, "A1 theta bin 0, energy bin 0");
  CheckClose(Factor(0.5f, 25.f), 1.2f, "A1 theta bin 0, energy bin 1");
  CheckClose(Factor(0.5f, 75.f), 1.3f, "A1 theta bin 0, energy bin 2");
  CheckClose(Factor(1.5f, 5.f), 2.1f, "A1 theta bin 1, energy bin 0");
  CheckClose(Factor(1.5f, 25.f), 2.2f, "A1 theta bin 1, energy bin 1");
  CheckClose(Factor(1.5f, 75.f), 2.3f, "A2 theta bin 1, energy bin 2 (catches transposed index)");

  CheckClose(Factor(0.5f, 10.f), 1.2f, "A3 energy exactly on an interior edge is in the upper bin");
  CheckClose(Factor(0.5f, 9.999f), 1.1f, "A4 energy just below an interior edge is in the lower bin");
  CheckClose(Factor(0.f, 25.f), 1.2f, "A5 value exactly on front() is in bin 0");

  // A6: theta cannot be placed exactly on a bin edge through this API. The caller
  // supplies a direction, and the lookup recovers theta as acos(cos(theta)); for
  // theta = 1.0f that round trip returns 0.99999994f, one ULP low, so the value
  // lands in the LOWER bin. This is a characterisation, not a requirement -- it is
  // the reason the boundary behaviour in theta is untestable until the lookup
  // exposes a GetCorrection(theta, energy) taking theta directly.
  CheckClose(Factor(1.0f, 25.f), 1.2f,
             "A6 direction built at exactly a theta edge lands in the LOWER bin (acos/cos round trip)");
}

// B: out-of-range behaviour.
void CaseOutOfRange(const Pandora& pandora) {
  RegisterAsymmetricTable(pandora);

  CheckClose(Factor(0.5f, 150.f), 1.0f, "B1 energy above back() -> no correction");
  CheckClose(Factor(0.5f, 100.f), 1.0f, "B2 energy exactly at back() -> no correction (>= not >)");
  CheckClose(Factor(0.5f, 99.999f), 1.3f, "B3 energy just below back() -> still corrected");
  CheckClose(Factor(2.5f, 25.f), 1.0f, "B4 theta above back() -> no correction");
  CheckClose(Factor(2.0f, 25.f), 1.0f, "B5 theta exactly at back() -> no correction");

  // The cliff, stated as a ratio. This pins CURRENT behaviour; clamping to the edge
  // bins would make the ratio 1.0 and this assertion is the one to flip.
  const float below(LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, DirAtTheta(0.5f), 99.999f));
  const float above(LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, DirAtTheta(0.5f), 100.001f));
  std::cout << "  note  B6 discontinuity at the top energy edge: " << below << " -> " << above << " GeV ("
            << (above / below - 1.f) * 100.f << "% jump for 2 MeV of raw energy)" << std::endl;
  Check(above < below, "B6 corrected energy DECREASES across the top edge (non-monotonic, factor > 1)");
}

// C: direction handling.
void CaseDirection(const Pandora& pandora) {
  const FloatVector thetaEdges{0.f, 1.f, 2.f, 3.2f};
  const FloatVector energyEdges{0.f, 100.f};
  const FloatVector factors{1.1f, 1.2f, 1.3f};
  if (STATUS_CODE_SUCCESS !=
      LCContent::RegisterNonLinearityEnergyCorrection(pandora, "Dir", HADRONIC, thetaEdges, energyEdges, factors))
    throw StatusCodeException(STATUS_CODE_FAILURE);

  CheckClose(Factor(0.f, 50.f), 1.1f, "C1 +z direction -> theta 0");
  CheckClose(Factor(static_cast<float>(M_PI) / 2.f, 50.f), 1.2f, "C2 transverse direction -> theta pi/2");
  CheckClose(Factor(static_cast<float>(M_PI), 50.f), 1.3f, "C3 -z direction -> theta pi");

  const float unit(
      LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, CartesianVector(0.f, 0.f, 1.f), 50.f));
  const float scaled(
      LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, CartesianVector(0.f, 0.f, 500.f), 50.f));
  CheckClose(scaled, unit, "C4 direction magnitude does not affect the result");

  const float zero(
      LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, CartesianVector(0.f, 0.f, 0.f), 50.f));
  CheckClose(zero, 50.f, "C5 zero-magnitude direction returns the input energy");
}

// D: rejection of malformed tables.
void CaseValidation(const Pandora& pandora) {
  ExpectRejected(pandora, "D1", {0.f, 1.f, 2.f}, {0.f, 10.f, 50.f}, {1.f, 1.f, 1.f}, "D1 wrong number of factors");
  ExpectRejected(pandora, "D2", {0.f, 2.f, 1.f}, {0.f, 10.f}, {1.f, 1.f}, "D2 non-monotonic theta edges");
  ExpectRejected(pandora, "D3", {0.f, 1.f, 1.f}, {0.f, 10.f}, {1.f, 1.f}, "D3 duplicate adjacent theta edges");
  ExpectRejected(pandora, "D4", {0.f}, {0.f, 10.f}, {1.f}, "D4 fewer than two theta edges");
  ExpectRejected(pandora, "D5", {0.f, 1.f}, {0.f}, {1.f}, "D5 fewer than two energy edges");
}

// E: the advertised no-op guarantee.
void CaseNoTable(const Pandora&) {
  CheckClose(LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, DirAtTheta(1.f), 42.f), 42.f,
             "E1 no registered table -> energy returned unchanged");
}

void CaseTypeIsolation(const Pandora& pandora) {
  const FloatVector thetaEdges{0.f, 2.f};
  const FloatVector energyEdges{0.f, 100.f};
  const FloatVector factors{3.f};
  if (STATUS_CODE_SUCCESS != LCContent::RegisterNonLinearityEnergyCorrection(pandora, "EmOnly", ELECTROMAGNETIC,
                                                                             thetaEdges, energyEdges, factors))
    throw StatusCodeException(STATUS_CODE_FAILURE);

  CheckClose(LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(HADRONIC, DirAtTheta(1.f), 50.f), 50.f,
             "E2 an ELECTROMAGNETIC table is not applied to a HADRONIC query");
  CheckClose(LCEnergyCorrectionPlugins::GetThetaEnergyCorrectedEnergy(ELECTROMAGNETIC, DirAtTheta(1.f), 50.f), 150.f,
             "E2 the ELECTROMAGNETIC table is applied to an ELECTROMAGNETIC query");
}

// F: requirements on how the registry resolves a table. Both cases assert the
// behaviour the lookup should have; neither holds today, because
// GetThetaEnergyCorrectedEnergy is keyed on the correction type alone and has no
// way to be told which table the caller wants.
void CaseNameAmbiguity(const Pandora& pandora) {
  const FloatVector thetaEdges{0.f, 2.f};
  const FloatVector energyEdges{0.f, 100.f};
  LCContent::RegisterNonLinearityEnergyCorrection(pandora, "AAA", HADRONIC, thetaEdges, energyEdges, {1.5f});
  LCContent::RegisterNonLinearityEnergyCorrection(pandora, "ZZZ", HADRONIC, thetaEdges, energyEdges, {2.5f});

  // With two tables registered for one correction type and no name supplied at the
  // call, the request is ambiguous. Applying either factor is a silent guess whose
  // outcome depends on how std::map happens to order the two name strings, so the
  // lookup must decline and leave the energy alone.
  //
  // Once the lookup takes the plugin name, replace this with a pair of assertions
  // that each name resolves to its own factor.
  const float f(Factor(1.f, 50.f));
  std::cout << "  note  F1 two HADRONIC tables registered (AAA=1.5, ZZZ=2.5); lookup returned " << f << std::endl;
  CheckClose(f, 1.0f, "F1 an ambiguous lookup must not silently apply an arbitrarily chosen table");
}

void CaseMultiInstance(const Pandora&) {
  const Pandora pandoraA, pandoraB;
  const FloatVector thetaEdges{0.f, 2.f};
  const FloatVector energyEdges{0.f, 100.f};

  LCContent::RegisterNonLinearityEnergyCorrection(pandoraA, "Shared", HADRONIC, thetaEdges, energyEdges, {1.5f});
  const float afterA(Factor(1.f, 50.f));
  CheckClose(afterA, 1.5f, "F2 a table registered on the first instance is visible");

  // Registering on a second, independent Pandora instance must not disturb what the
  // first instance had. Each instance carries its own calibration.
  LCContent::RegisterNonLinearityEnergyCorrection(pandoraB, "Shared", HADRONIC, thetaEdges, energyEdges, {2.5f});
  const float afterB(Factor(1.f, 50.f));
  std::cout << "  note  F2 instance A registered 1.5, instance B registered 2.5; lookup returned " << afterB
            << std::endl;
  CheckClose(afterB, 1.5f, "F2 registering on a second instance must not overwrite the first instance's table");
}

} // namespace

int main(int argc, char* argv[]) {
  const std::string testCase(argc > 1 ? argv[1] : "");

  try {
    const Pandora pandora;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(pandora));

    std::cout << "case: " << testCase << std::endl;

    if ("lookup" == testCase)
      CaseLookup(pandora);
    else if ("out-of-range" == testCase)
      CaseOutOfRange(pandora);
    else if ("direction" == testCase)
      CaseDirection(pandora);
    else if ("validation" == testCase)
      CaseValidation(pandora);
    else if ("no-table" == testCase)
      CaseNoTable(pandora);
    else if ("type-isolation" == testCase)
      CaseTypeIsolation(pandora);
    else if ("name-ambiguity" == testCase)
      CaseNameAmbiguity(pandora);
    else if ("multi-instance" == testCase)
      CaseMultiInstance(pandora);
    else {
      std::cerr << "unknown test case '" << testCase << "'" << std::endl;
      return 2;
    }
  } catch (const StatusCodeException& e) {
    std::cerr << "  ERROR unexpected pandora exception: " << e.ToString() << std::endl;
    return 1;
  }

  std::cout << (g_failures ? "RESULT: FAILED" : "RESULT: passed") << " (" << g_failures << " failure(s))" << std::endl;
  return g_failures ? 1 : 0;
}
