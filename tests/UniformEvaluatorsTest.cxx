/*!
 * \file   UniformEvaluatorsTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <memory>
#include <limits>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/UniformEvaluator.hxx"
#include "MGIS/Function/Tensors.hxx"

struct UniformEvaluatorsTest final : public tfel::tests::TestCase {
  UniformEvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function", "UniformEvaluatorsTests") {
  }  // end of UniformEvaluatorsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) constexpr {
      return r > 0 ? r : -r;
    };
    constexpr auto eps = real{1e-14};
    constexpr auto space = BasicLinearSpace{3};
    constexpr auto e = UniformEvaluator(space, 12.);
    TFEL_TESTS_STATIC_ASSERT(local_abs(e(1) - 12) < eps);
    constexpr real values[3] = {12, -2, 6};
    constexpr auto e2 = UniformEvaluator(space, values);
    TFEL_TESTS_STATIC_ASSERT(local_abs(e2(1)[0] - 12) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(e2(1)[1] + 2) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(e2(1)[2] - 6) < eps);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) constexpr {
      return r > 0 ? r : -r;
    };
    constexpr auto eps = real{1e-14};
    constexpr auto space = BasicLinearSpace{3};
    constexpr real values[6] = {1, 1, 1, 0, 0, 0};
    constexpr auto id = UniformEvaluator(space, values) | as_stensor<3>;
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[0] - 1) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[1] - 1) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[2] - 1) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[3]) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[4]) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(id(1)[5]) < eps);
  }
};

TFEL_TESTS_GENERATE_PROXY(UniformEvaluatorsTest, "UniformEvaluatorsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("UniformEvaluatorsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
