/*!
 * \file   MechanicalEvaluatorsTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <memory>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/FixedSizeEvaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators.hxx"

struct MechanicalEvaluatorsTest final : public tfel::tests::TestCase {
  MechanicalEvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function", "MechanicalEvaluatorsTests") {
  }  // end of MechanicalEvaluatorsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(1);
    auto values = std::vector<real>{1, 2, 3};
    const auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 3);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    const auto s = FixedSizeEvaluator<BasicLinearSpace, 3>(f);
    const auto e = vmis<1>(s);
    TFEL_TESTS_ASSERT(std::abs(e(0) - std::sqrt(3)) < eps);
    const auto e2 = hydrostatic_stress<1>(s);
    TFEL_TESTS_ASSERT(std::abs(e2(0) - 2) < eps);
  }
};

TFEL_TESTS_GENERATE_PROXY(MechanicalEvaluatorsTest, "MechanicalEvaluatorsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("MechanicalEvaluatorsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
