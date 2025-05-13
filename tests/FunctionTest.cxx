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
#include "MGIS/Function/FixedSizeView.hxx"

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
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    auto space = std::make_shared<BasicLinearSpace>(0);
    auto values = std::vector<real>{};
    Context ctx;
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            ctx, space, values, 2, 3));
    auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 2, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3};
    const auto ok = ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
        ctx, space, values, 1);
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
            ctx, space, values2, 1));
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
            ctx, space, values3, 1));
    TFEL_TESTS_CHECK_THROW(
        ImmutableFunctionView<BasicLinearSpace>(space, values3, 1),
        std::runtime_error);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(ctx, space,
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
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(3);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    const auto ok =
        ImmutableFunctionView<BasicLinearSpace,
                              {.data_size = 2}>::checkPreconditions(ctx, space,
                                                               values, 2);
    TFEL_TESTS_ASSERT(ok);
    auto f =
        ImmutableFunctionView<BasicLinearSpace, {.data_size = 2}>(space, values, 2);
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
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<BasicLinearSpace>::checkPreconditions(
            ctx, space, values, 2, 3));
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
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        (ImmutableFunctionView<
            BasicLinearSpace,
            {.data_size = 2, .data_stride = 3}>::checkPreconditions(ctx, space,
                                                                    values)));
    auto f = ImmutableFunctionView<BasicLinearSpace, {.data_size = 2, .data_stride = 3}>(
        space, values);
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_STATIC_ASSERT(f.getDataStride() == 3);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.getValues(1)[1] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f(1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f(1)[1] - 5) < eps);
  }
  void test7() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearQuadratureSpace<3>>(2);
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        ImmutableFunctionView<
            BasicLinearQuadratureSpace<3>>::checkPreconditions(ctx, space,
                                                               values, 1, 1));
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
};

TFEL_TESTS_GENERATE_PROXY(ImmutableFunctionTest, "ImmutableFunctionTest");

struct FunctionTest final : public tfel::tests::TestCase {
  FunctionTest()
      : tfel::tests::TestCase("MGIS/Function", "FunctionTest") {
  }  // end of FunctionTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    this->test6();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f = Function<BasicLinearSpace>(space, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    f.getValues(0)[0] = 1;
    f.getValues(0)[1] = 2;
    f.getValues(1)[0] = 3;
    f.getValues(1)[1] = 4;
    const auto values = f.getValues();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[3] - 4) < eps);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f = Function<BasicLinearSpace>(space, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    f(0)[0] = 1;
    f(0)[1] = 2;
    f(1)[0] = 3;
    f(1)[1] = 4;
    const auto values = f.getValues();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[3] - 4) < eps);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f2 = [&space] {
      auto f = Function<BasicLinearSpace>(space, 2);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    TFEL_TESTS_ASSERT(space.get() == &(f2.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f2.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f2.getDataStride(), 2);
    const auto values2 = f2.getValues();
    TFEL_TESTS_ASSERT(std::abs(values2[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[3] - 4) < eps);
    // move constructor
    auto f3 = std::move(f2);
    TFEL_TESTS_ASSERT(space.get() == &(f3.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f3.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f3.getDataStride(), 2);
    const auto values3 = f3.getValues();
    TFEL_TESTS_ASSERT(std::abs(values3[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[3] - 4) < eps);
    // copy constructor
    auto f4(f3);
    TFEL_TESTS_ASSERT(space.get() == &(f4.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f4.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f4.getDataStride(), 2);
    const auto values4 = f4.getValues();
    TFEL_TESTS_ASSERT(std::abs(values4[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[3] - 4) < eps);
    // view
    auto f5 = view(f4);
    TFEL_TESTS_ASSERT(space.get() == &(f5.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f5.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f5.getDataStride(), 2);
    const auto values5 = f5.getValues();
    TFEL_TESTS_ASSERT(std::abs(values5[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values5[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values5[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values5[3] - 4) < eps);
    //
  }
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f2 = [&space] {
      auto f = Function<BasicLinearSpace>(space, 2);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    auto f3 = Function<BasicLinearSpace>(space, 2);
    static_assert(std::same_as<decltype(view(f2)),
                               FunctionView<BasicLinearSpace, {}, false>>);
    const auto ok = assign(ctx, f3, view(f2));
    TFEL_TESTS_ASSERT(ok);
    const auto values3 = f3.getValues();
    TFEL_TESTS_ASSERT(std::abs(values3[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[3] - 4) < eps);
    auto f4 = Function<BasicLinearSpace>(space, 2);
    const auto ok2 = view(f2) | f4;
    TFEL_TESTS_ASSERT(ok2);
    const auto values4 = f4.getValues();
    TFEL_TESTS_ASSERT(std::abs(values4[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[3] - 4) < eps);
  }
  void test5() {
    // check that view<N> work with a dynamic sized function
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f2 = [&space] {
      auto f = Function<BasicLinearSpace>(space, 2);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    auto f3 = view<2>(f2);
    TFEL_TESTS_ASSERT(f3.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(f3(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(1)[1] - 4) < eps);
    // changing f2 value using the view
    f3(0)[0] = 4;
    f3(0)[1] = 5;
    f3(1)[0] = 6;
    f3(1)[1] = 9;
    TFEL_TESTS_ASSERT(std::abs(f2(0)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f2(0)[1] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f2(1)[0] - 6) < eps);
    TFEL_TESTS_ASSERT(std::abs(f2(1)[1] - 9) < eps);
    // scalar case
    auto f4 = [&space] {
      auto f = Function<BasicLinearSpace>(space, 1);
      f(0)[0] = 1;
      f(1)[0] = 2;
      return f;
    }();
    auto f5 = view<1>(f4);
    TFEL_TESTS_ASSERT(f5.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(f5(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f5(1) - 2) < eps);
    // changing f4 value using the view
    f5(0) = 4;
    f5(1) = 5;
    TFEL_TESTS_ASSERT(std::abs(f4(0)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f4(1)[0] - 5) < eps);
  }
  void test6() {
    // check that view<N> work with a fixed size function
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto f2 = [&space] {
      auto f = Function<BasicLinearSpace, 2>(space);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    const auto f3 = view<2>(f2);
    TFEL_TESTS_ASSERT(std::abs(f3(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f3(1)[1] - 4) < eps);
  }
};

TFEL_TESTS_GENERATE_PROXY(FunctionTest, "FunctionTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FunctionTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
