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
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3};
    const auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(1) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(2) - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[0] - 3) < eps);
    auto values2 = std::vector<real>{1, 2, 3, 4};
    const auto f2 = ImmutableFunctionView<BasicLinearSpace>(space, values2, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(1) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(2) - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[0] - 3) < eps);
    auto values3 = std::vector<real>{1, 2};
    TFEL_TESTS_CHECK_THROW(
        ImmutableFunctionView<BasicLinearSpace>(space, values3, 1),
        std::runtime_error);
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
