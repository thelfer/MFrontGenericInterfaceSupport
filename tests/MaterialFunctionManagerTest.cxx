/*!
 * \file   MaterialFunctionManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   03/06/2025
 * \copyright (C) Copyright Thomas Helfer 2025.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialFunctionManager.hxx"

static const char* library = nullptr;

struct MaterialFunctionManagerTest final : public tfel::tests::TestCase {
  MaterialFunctionManagerTest()
      : tfel::tests::TestCase("MGIS/Function", "MaterialFunctionManagerTest") {
  }  // end of MechanicalEvaluatorsTest
  tfel::tests::TestResult execute() override {
    this->test1<true>();
    this->test1<false>();
    return this->result;
  }

 private:
  template <bool mutable_functions>
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    using namespace mgis::behaviour;
    constexpr auto eps = real{1e-15};
    using FunctionManager =
        std::conditional_t<mutable_functions,
                           MaterialFunctionManager<BasicLinearSpace>,
                           const MaterialFunctionManager<BasicLinearSpace>>;
    const auto b = load(library, "Norton", Hypothesis::TRIDIMENSIONAL);
    auto qspace = std::make_shared<BasicLinearSpace>(3);
    FunctionManager m = MaterialFunctionManager<BasicLinearSpace>{qspace, b};
    auto ctx = Context{};
    /* gradients */
    auto oeto = getGradient(ctx, m, "Strain");
    auto oeto2 = getGradient(ctx, m, "Strain2");
    auto osig = getThermodynamicForce(ctx, m, "Stress");
    auto osig2 = getThermodynamicForce(ctx, m, "Stress2");
    TFEL_TESTS_ASSERT(isValid(oeto));
    TFEL_TESTS_CHECK_EQUAL(&(*(oeto->getSpace())), qspace.get());
    TFEL_TESTS_CHECK_EQUAL(oeto->getNumberOfComponents(), 6);
    TFEL_TESTS_ASSERT(isInvalid(oeto2));
    TFEL_TESTS_ASSERT(isValid(osig));
    TFEL_TESTS_CHECK_EQUAL(&(*(osig->getSpace())), qspace.get());
    TFEL_TESTS_CHECK_EQUAL(osig->getNumberOfComponents(), 6);
    TFEL_TESTS_ASSERT(isInvalid(osig2));
    /* internal state variables */
    auto oeel = getInternalStateVariable(ctx, m, "ElasticStrain");
    TFEL_TESTS_ASSERT(isValid(oeel));
    TFEL_TESTS_CHECK_EQUAL(&(*(oeel->getSpace())), qspace.get());
    auto op = getInternalStateVariable(ctx, m, "EquivalentViscoplasticStrain");
    TFEL_TESTS_ASSERT(isValid(op));
    TFEL_TESTS_CHECK_EQUAL(&(*(op->getSpace())), qspace.get());
    if constexpr (mutable_functions) {
      auto value = [](const size_type i) -> real { return (2 * i + 1) * 1e-2; };
      for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
        (*op)(i)[0] = value(i);
      }
      const auto& isvs = m.s1.internal_state_variables;
      for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
        TFEL_TESTS_ASSERT(std::abs(isvs[i * 7 + 6] - value(i)) < eps);
      }
    }
    // inexistant internal state variabe variable
    auto op2 = getInternalStateVariable(ctx, m, "EquivalentPlasticStrain");
    TFEL_TESTS_ASSERT(isInvalid(op2));
    // access with known size
    auto op3 =
        getInternalStateVariable<1>(ctx, m, "EquivalentViscoplasticStrain");
    TFEL_TESTS_ASSERT(isValid(op3));
    TFEL_TESTS_CHECK_EQUAL(&(*(op3->getSpace())), qspace.get());
    if constexpr (mutable_functions) {
      auto value = [](const size_type i) -> real {
        return std::exp(2 * i) * 1e-2;
      };
      for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
        (*op3)(i) = value(i);
      }
      const auto& isvs = m.s1.internal_state_variables;
      for (size_type i = 0; i != getSpaceSize(*qspace); ++i) {
        TFEL_TESTS_ASSERT(std::abs(isvs[i * 7 + 6] - value(i)) < eps);
      }
    }
    // invalid size
    auto op4 =
        getInternalStateVariable<4>(ctx, m, "EquivalentViscoplasticStrain");
    TFEL_TESTS_ASSERT(isInvalid(op4));
  }
};

TFEL_TESTS_GENERATE_PROXY(MaterialFunctionManagerTest,
                          "MaterialFunctionManagerTest");

int main(const int argc, const char* const* argv) {
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  library = argv[1];
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("MaterialFunctionManagerTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
