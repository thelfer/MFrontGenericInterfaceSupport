/*!
 * \file   MechancialEvaluatorsTest.cxx
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
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/FixedSizeEvaluator.hxx"
#include "MGIS/Function/Tensors.hxx"

struct EvaluatorsTest final : public tfel::tests::TestCase {
  EvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function", "EvaluatorsTests") {
  }  // end of EvaluatorsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3};
    const auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    const auto e = FixedSizeEvaluator<BasicLinearSpace, 1>(f);
    TFEL_TESTS_ASSERT(e.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(e(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(e(1) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(e(2) - 3) < eps);
    const auto e2 = view<1>(f) | transform([](const real x) { return x * x; });
    TFEL_TESTS_ASSERT(e2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(e2(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(e2(1) - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(e2(2) - 9) < eps);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    struct InvalidEvaluator {
      double getSpace() const;
    };
    TFEL_TESTS_STATIC_ASSERT(!EvaluatorConcept<double>);
    TFEL_TESTS_STATIC_ASSERT(!EvaluatorConcept<Test>);
    TFEL_TESTS_STATIC_ASSERT(
        EvaluatorConcept<ImmutableFunctionView<BasicLinearSpace>>);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(1);
    auto values = std::vector<real>{1, 2, 3};
    const auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 3);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    const auto e = view<3>(f);
    TFEL_TESTS_ASSERT(e.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(e(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(e(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(e(0)[2] - 3) < eps);
    const auto e2 = e | as_stensor<1>;
    TFEL_TESTS_ASSERT(e2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(tfel::math::trace(e2(0)) - 6) < eps);
    const auto e3 = e | as_stensor<1> | trace;
    TFEL_TESTS_ASSERT(e3.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(e3(0) - 6) < eps);
  }
};

TFEL_TESTS_GENERATE_PROXY(EvaluatorsTest, "EvaluatorsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("EvaluatorsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
