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
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"

struct ImmutableFunctionTest final : public tfel::tests::TestCase {
  ImmutableFunctionTest()
      : tfel::tests::TestCase("MGIS/Function", "ImmutableFunctionTests") {
  }  // end of ImmutableFunctionTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    this->test6();
    this->test7();
    this->test8();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    auto space = std::make_shared<BasicLinearSpace>(0);
    auto values = std::vector<real>{};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            space, values, 2, 3));
    auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 2, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3};
    const auto ok = ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
        space, values, 1);
    TFEL_TESTS_ASSERT(ok);
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
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            space, values2, 1));
    const auto f2 = ImmutableFunctionView<BasicLinearSpace>(space, values2, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(1) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(2) - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[0] - 3) < eps);
    // precondition is not met
    auto values3 = std::vector<real>{1, 2};
    TFEL_TESTS_ASSERT(
        !ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            space, values3, 1));
    TFEL_TESTS_CHECK_THROW(
        ImmutableFunctionView<BasicLinearSpace>(space, values3, 1),
        std::runtime_error);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(space,
                                                                    values, 2));
    auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(1) - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(2) - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[1] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[0] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[1] - 6) < eps);
  }
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    const auto ok =
        ImmutableFunctionView<BasicLinearSpace,
                              {.size = 2}>::checkPreconditions(space, values,
                                                               2);
    TFEL_TESTS_ASSERT(ok);
    auto f =
        ImmutableFunctionView<BasicLinearSpace, {.size = 2}>(space, values, 2);
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[1] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[0] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(2)[1] - 6) < eps);
  }
  void test5() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            space, values, 2, 3));
    auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 2, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValue(1) - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[1] - 5) < eps);
  }
  void test6() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        (ImmutableFunctionView<BasicLinearSpace, {.size = 2, .stride = 3}>::
             checkPreconditions(space, values)));
    auto f = ImmutableFunctionView<BasicLinearSpace, {.size = 2, .stride = 3}>(
        space, values);
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_STATIC_ASSERT(f.getDataStride() == 3);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[1] - 5) < eps);
  }
  void test7() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearQuadratureSpace<3>>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(ImmutableFunctionView<
                      BasicLinearQuadratureSpace<3>>::checkPreconditions(space,
                                                                         values,
                                                                         1, 1));
    auto f = ImmutableFunctionView<BasicLinearQuadratureSpace<3>>(space, values,
                                                                  1, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    for (size_type i = 0; i != space->size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(f.getValue(i) - (i + 1)) < eps);
      TFEL_TESTS_ASSERT(std::abs(f.getValues(i)[0] - (i + 1)) < eps);
    }
    for (size_type e = 0; e != space->getNumberOfCells(); ++e) {
      for (size_type i = 0; i != space->getNumberOfQuadraturePoints(e); ++i) {
        const auto v = e * 3 + i + 1;
        TFEL_TESTS_ASSERT(std::abs(f.getValue(e, i) - v) < eps);
        TFEL_TESTS_ASSERT(std::abs(f.getValues(e, i)[0] - v) < eps);
      }
    }
  }
  void test8() {
    using namespace mgis;
    using namespace mgis::function;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f = Function<BasicLinearSpace>(space, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
  }
};

TFEL_TESTS_GENERATE_PROXY(ImmutableFunctionTest, "ImmutableFunctionTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FunctionTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
