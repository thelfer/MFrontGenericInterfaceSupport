/*!
 * \file   FunctionTest.cxx
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
#include "MGIS/Function/SharedSpace.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/FixedSizeView.hxx"

struct ReduceTest final : public tfel::tests::TestCase {
  ReduceTest()
      : tfel::tests::TestCase("MGIS/Function", "ReduceTests") {
  }  // end of ReduceTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(3)};
    auto values = std::vector<real>{1, -1, 2};
    Context ctx;
    auto f = FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    const auto max_value = [](const real a, const real b) {
      return std::max(a, b);
    };
    const auto oresult = scalar_reduce(ctx, f, max_value, 0);
    TFEL_TESTS_ASSERT(isValid(oresult));
    std::cout << "value: " << *oresult << '\n';
    TFEL_TESTS_ASSERT(std::abs(*oresult - 2) < 1e-14);
  }
};

TFEL_TESTS_GENERATE_PROXY(ReduceTest, "ReduceTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ReduceTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
