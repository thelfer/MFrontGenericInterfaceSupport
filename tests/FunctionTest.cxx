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
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(0)};
    auto values = std::vector<real>{};
    Context ctx;
    TFEL_TESTS_ASSERT(
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values, 2, 3));
    auto f =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values, 2, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(3)};
    auto values = std::vector<real>{1, 2, 3};
    const auto ok =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values, 1);
    TFEL_TESTS_ASSERT(ok);
    const auto f =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[0] - 3) < eps);
    auto values2 = std::vector<real>{1, 2, 3, 4};
    TFEL_TESTS_ASSERT(
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values2, 1));
    const auto f2 =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values2, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[0] - 3) < eps);
    // precondition is not met
    auto values3 = std::vector<real>{1, 2};
#ifndef MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
    TFEL_TESTS_ASSERT(
        !FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values3, 1));
#endif MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
#ifdef MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION
    TFEL_TESTS_CHECK_THROW(
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values3, 1),
        std::runtime_error);
#endif /* MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION */
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(3)};
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values, 2));
    auto f = FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[1] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[0] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[1] - 6) < eps);
  }
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(3)};
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    const auto ok =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>,
                          {.data_size = 2}>::checkPreconditions(ctx, space,
                                                                values, 2);
    TFEL_TESTS_ASSERT(ok);
    auto f = FunctionEvaluator<SharedSpace<BasicLinearSpace>, {.data_size = 2}>(
        space, values, 2);
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[1] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[0] - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 2)[1] - 6) < eps);
  }
  void test5() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>::checkPreconditions(
            ctx, space, values, 2, 3));
    auto f =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>>(space, values, 2, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[1] - 5) < eps);
  }
  void test6() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        (FunctionEvaluator<SharedSpace<BasicLinearSpace>,
                           {.data_size = 2,
                            .data_stride = 3}>::checkPreconditions(ctx, space,
                                                                   values)));
    auto f =
        FunctionEvaluator<SharedSpace<BasicLinearSpace>,
                          {.data_size = 2, .data_stride = 3}>(space, values);
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_STATIC_ASSERT(f.getDataStride() == 3);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, 1)[1] - 5) < eps);
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
    auto space =
        SharedSpace{std::make_shared<BasicLinearQuadratureSpace<3>>(2)};
    auto values = std::vector<real>{1, 2, 3, 4, 5, 6};
    TFEL_TESTS_ASSERT(
        FunctionEvaluator<SharedSpace<BasicLinearQuadratureSpace<3>>>::
            checkPreconditions(ctx, space, values, 1, 1));
    auto f = FunctionEvaluator<SharedSpace<BasicLinearQuadratureSpace<3>>>(
        space, values, 1, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    for (size_type i = 0; i != getSpaceSize(space); ++i) {
      TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, i)[0] - (i + 1)) < eps);
    }
    for (size_type e = 0; e != getNumberOfCells(space); ++e) {
      for (size_type i = 0; i != getNumberOfQuadraturePoints(space, e); ++i) {
        const auto v = e * 3 + i + 1;
        TFEL_TESTS_ASSERT(std::abs(f.data(unsafe, e, i)[0] - v) < eps);
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
    this->test7();
    this->test8();
    this->test9();
    this->test10();
    this->test11();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f = Function<SharedSpace<BasicLinearSpace>>(space, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    f.data(unsafe, 0)[0] = 1;
    f.data(unsafe, 0)[1] = 2;
    f.data(unsafe, 1)[0] = 3;
    f.data(unsafe, 1)[1] = 4;
    const auto values = f.data();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[3] - 4) < eps);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f = Function<SharedSpace<BasicLinearSpace>>(space, 2);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 2);
    f(0)[0] = 1;
    f(0)[1] = 2;
    f(1)[0] = 3;
    f(1)[1] = 4;
    const auto values = f.data();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[3] - 4) < eps);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f2 = [&space] {
      auto f = Function<SharedSpace<BasicLinearSpace>>(space, 2);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    TFEL_TESTS_ASSERT(areEquivalent(space, f2.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f2.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f2.getDataStride(), 2);
    const auto values2 = f2.data();
    TFEL_TESTS_ASSERT(std::abs(values2[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values2[3] - 4) < eps);
    // move constructor
    auto f3 = std::move(f2);
    TFEL_TESTS_ASSERT(areEquivalent(space, f3.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f3.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f3.getDataStride(), 2);
    const auto values3 = f3.data();
    TFEL_TESTS_ASSERT(std::abs(values3[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[3] - 4) < eps);
    // copy constructor
    auto f4(f3);
    TFEL_TESTS_ASSERT(areEquivalent(space, f4.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f4.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f4.getDataStride(), 2);
    const auto values4 = f4.data();
    TFEL_TESTS_ASSERT(std::abs(values4[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values4[3] - 4) < eps);
    // view
    auto f5 = view(f4);
    TFEL_TESTS_ASSERT(areEquivalent(space, f5.getSpace()));
    TFEL_TESTS_CHECK_EQUAL(f5.getNumberOfComponents(), 2);
    TFEL_TESTS_CHECK_EQUAL(f5.getDataStride(), 2);
    const auto values5 = f5.data();
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
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f2 = [&space] {
      auto f = Function<SharedSpace<BasicLinearSpace>>(space, 2);
      f(0)[0] = 1;
      f(0)[1] = 2;
      f(1)[0] = 3;
      f(1)[1] = 4;
      return f;
    }();
    auto f3 = Function<SharedSpace<BasicLinearSpace>>(space, 2);
    static_assert(
        std::same_as<decltype(view(f2)),
                     FunctionView<SharedSpace<BasicLinearSpace>, {}, false>>);
    const auto ok = assign(ctx, f3, view(f2));
    TFEL_TESTS_ASSERT(ok);
    const auto values3 = f3.data();
    TFEL_TESTS_ASSERT(std::abs(values3[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[2] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values3[3] - 4) < eps);
    auto f4 = Function<SharedSpace<BasicLinearSpace>>(space, 2);
    const auto ok2 = view(f2) | f4;
    TFEL_TESTS_ASSERT(ok2);
    const auto values4 = f4.data();
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
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f2 = [&space] {
      auto f = Function<SharedSpace<BasicLinearSpace>>(space, 2);
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
      auto f = Function<SharedSpace<BasicLinearSpace>>(space, 1);
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
    //
    auto f6 = f4 | as_scalar;
    TFEL_TESTS_ASSERT(f6.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(f6(0) - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(f6(1) - 5) < eps);
    // changing f4 value using the view
    f6(0) = 2;
    f6(1) = -3;
    TFEL_TESTS_ASSERT(std::abs(f4(0)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(f4(1)[0] + 3) < eps);
  }
  void test6() {
    // check that view<N> work with a fixed size function
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f2 = [&space] {
      auto f = Function<SharedSpace<BasicLinearSpace>, 2>(space);
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
  void test7() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = SharedSpace{std::make_shared<BasicLinearSpace>(2)};
    auto f = Function<SharedSpace<BasicLinearSpace>, 1>(space);
    TFEL_TESTS_STATIC_ASSERT(FunctionConcept<decltype(f)>);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    TFEL_TESTS_STATIC_ASSERT(number_of_components<decltype(f)> == 1);
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<function_result<decltype(f)>, real&>));
    f(0) = 12;
    f(1) = 13;
    const auto& values = f.data();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 12) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 13) < eps);
    const auto& f2 = f;
    TFEL_TESTS_ASSERT(std::abs(f2(0) - 12) < eps);
    TFEL_TESTS_ASSERT(std::abs(f2(1) - 13) < eps);
  }
  void test8() {
    using namespace mgis;
    using namespace mgis::function;
    auto check_value = [](const real& a, const real b) constexpr->bool {
      constexpr auto eps = real{1e-12};
      auto local_abs = [](const real r) { return r > 0 ? r : -r; };
      return local_abs(a - b) < eps;
    };
    static constexpr auto space = BasicLinearSpace{2};
    static constexpr auto values = std::array<real, 4>{4, 2, 3, 1};
    constexpr auto f = FunctionEvaluator<BasicLinearSpace>{space, values, 2};
    TFEL_TESTS_STATIC_ASSERT(f.getNumberOfComponents() == 2);
    TFEL_TESTS_STATIC_ASSERT(f.getDataStride() == 2);
    TFEL_TESTS_STATIC_ASSERT(check_value(f(0)[0], 4));
    TFEL_TESTS_STATIC_ASSERT(check_value(f(0)[1], 2));
    TFEL_TESTS_STATIC_ASSERT(check_value(f(1)[0], 3));
    TFEL_TESTS_STATIC_ASSERT(check_value(f(1)[1], 1));
  }
  void test9() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    auto check_value = [](const real& a, const real b) constexpr->bool {
      constexpr auto eps = real{1e-12};
      auto local_abs = [](const real r) { return r > 0 ? r : -r; };
      return local_abs(a - b) < eps;
    };
    constexpr auto sizes = []() constexpr {
      auto space = BasicLinearSpace{2};
      auto f = Function<BasicLinearSpace>{space, 3};
      return std::array<size_type, 2>{f.getNumberOfComponents(),
                                      f.getDataStride()};
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(sizes[0] == 3);
    TFEL_TESTS_STATIC_ASSERT(sizes[1] == 3);
    constexpr auto values = []() constexpr->std::array<real, 4> {
      auto space = BasicLinearSpace{2};
      auto f = Function<BasicLinearSpace>{space, 2};
      f(0)[0] = 5;
      f(0)[1] = 12;
      f(1)[0] = -2;
      f(1)[1] = 3;
      const auto data = f.data();
      std::array<real, 4> fvalues;
      std::copy(data.begin(), data.end(), fvalues.begin());
      return fvalues;
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(values[0], 5));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[1], 12));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[2], -2));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[3], 3));
    constexpr auto values2 = []() constexpr->std::array<real, 4> {
      auto space = BasicLinearSpace{2};
      auto f = Function<BasicLinearSpace>{space, 2};
      f.data(unsafe, 0)[0] = 5;
      f.data(unsafe, 0)[1] = -6;
      f.data(unsafe, 1)[0] = -2;
      f.data(unsafe, 1)[1] = 3;
      const auto data = f.data();
      std::array<real, 4> fvalues;
      std::copy(data.begin(), data.end(), fvalues.begin());
      return fvalues;
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(values2[0], 5));
    TFEL_TESTS_STATIC_ASSERT(check_value(values2[1], -6));
    TFEL_TESTS_STATIC_ASSERT(check_value(values2[2], -2));
    TFEL_TESTS_STATIC_ASSERT(check_value(values2[3], 3));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test10() {
    using namespace mgis;
    using namespace mgis::function;
#ifndef MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
    const auto s = BasicLinearSpace{2};
    auto values = std::vector<real>(5);
    auto check = [this, &s, &values](const size_type dsize, const auto dstride,
                                     std::string_view emsg) {
      auto ctx = Context{};
      const auto b = FunctionView<BasicLinearSpace>::checkPreconditions(
          ctx, s, values, dsize, dstride);
      TFEL_TESTS_ASSERT(!b);
      TFEL_TESTS_ASSERT(ctx.getRawErrorMessage() == emsg);
    };
    check(0, 2, "invalid number of components");
    check(2, 0, "invalid stride");
    check(3, 2, "the number of components is greater than the stride");
    check(3, 3, "insufficient external data size");
    auto check2 = [this, &s, &values](const size_type dsize,
                                      std::string_view emsg) {
      constexpr auto layout = FunctionDataLayoutDescription{
          .data_size = dynamic_extent, .data_stride = 3};
      auto ctx = Context{};
      const auto b = FunctionView<BasicLinearSpace, layout>::checkPreconditions(
          ctx, s, values, dsize);
      TFEL_TESTS_ASSERT(!b);
      TFEL_TESTS_ASSERT(ctx.getRawErrorMessage() == emsg);
    };
    check2(0, "invalid number of components");
    check2(4, "the number of components is greater than the stride");
    check2(3, "insufficient external data size");
    auto check3 = [this, &s, &values](const size_type dstride,
                                      std::string_view emsg) {
      constexpr auto layout = FunctionDataLayoutDescription{
          .data_size = 3, .data_stride = dynamic_extent};
      auto ctx = Context{};
      const auto b = FunctionView<BasicLinearSpace, layout>::checkPreconditions(
          ctx, s, values, dstride);
      TFEL_TESTS_ASSERT(!b);
      TFEL_TESTS_ASSERT(ctx.getRawErrorMessage() == emsg);
    };
    check3(0, "invalid stride");
    check3(2, "the number of components is greater than the stride");
    check3(3, "insufficient external data size");
    auto check4 = [this, &s, &values](std::string_view emsg) {
      constexpr auto layout =
          FunctionDataLayoutDescription{.data_size = 2, .data_stride = 4};
      auto ctx = Context{};
      const auto b = FunctionView<BasicLinearSpace, layout>::checkPreconditions(
          ctx, s, values);
      TFEL_TESTS_ASSERT(!b);
      TFEL_TESTS_ASSERT(ctx.getRawErrorMessage() == emsg);
    };
    check4("insufficient external data size");
    auto check5 = [this, &s](const size_type dsize, std::string_view emsg) {
      auto ctx = Context{};
      const auto b =
          Function<BasicLinearSpace>::checkPreconditions(ctx, s, dsize);
      TFEL_TESTS_ASSERT(!b);
      TFEL_TESTS_ASSERT(ctx.getRawErrorMessage() == emsg);
    };
    check5(0, "invalid number of components");
#endif /* MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING */
  }
  void test11() {
    using namespace mgis;
    using namespace mgis::function;
#ifdef MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION
    const auto s = BasicLinearSpace{2};
    auto values = std::vector<real>(5);
    auto check = [this, &s, &values](const size_type dsize, const auto dstride,
                                     std::string_view emsg) {
      auto exception_thrown = false;
      try {
        check_preconditions<FunctionView<BasicLinearSpace>>(s, values, dsize,
                                                            dstride);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
      exception_thrown = false;
      try {
        FunctionView<BasicLinearSpace> f(s, values, dsize, dstride);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
    };
    check(0, 2, "invalid number of components");
    check(2, 0, "invalid stride");
    check(3, 2, "the number of components is greater than the stride");
    check(3, 3, "insufficient external data size");
    auto check2 = [this, &s, &values](const size_type dsize,
                                      std::string_view emsg) {
      constexpr auto layout = FunctionDataLayoutDescription{
          .data_size = dynamic_extent, .data_stride = 3};
      auto exception_thrown = false;
      try {
        check_preconditions<FunctionView<BasicLinearSpace, layout>>(s, values,
                                                                    dsize);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
      exception_thrown = false;
      try {
        FunctionView<BasicLinearSpace, layout> f(s, values, dsize);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
    };
    check2(0, "invalid number of components");
    check2(4, "the number of components is greater than the stride");
    check2(3, "insufficient external data size");
    auto check3 = [this, &s, &values](const size_type dstride,
                                      std::string_view emsg) {
      constexpr auto layout = FunctionDataLayoutDescription{
          .data_size = 3, .data_stride = dynamic_extent};
      auto exception_thrown = false;
      try {
        check_preconditions<FunctionView<BasicLinearSpace, layout>>(s, values,
                                                                    dstride);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
      exception_thrown = false;
      try {
        FunctionView<BasicLinearSpace, layout> f(s, values, dstride);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
    };
    check3(0, "invalid stride");
    check3(2, "the number of components is greater than the stride");
    check3(3, "insufficient external data size");
    auto check4 = [this, &s, &values](std::string_view emsg) {
      constexpr auto layout =
          FunctionDataLayoutDescription{.data_size = 2, .data_stride = 4};
      auto exception_thrown = false;
      try {
        check_preconditions<FunctionView<BasicLinearSpace, layout>>(s, values);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
      exception_thrown = false;
      try {
        FunctionView<BasicLinearSpace, layout> f(s, values);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
    };
    check4("insufficient external data size");
    auto check5 = [this, &s, &values](const size_type dsize,
                                      std::string_view emsg) {
      auto exception_thrown = false;
      try {
        check_preconditions<Function<BasicLinearSpace>>(s, dsize);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
      exception_thrown = false;
      try {
        Function<BasicLinearSpace> f(s, dsize);
      } catch (std::exception& e) {
        TFEL_TESTS_ASSERT(e.what() == emsg);
        exception_thrown = true;
      }
      TFEL_TESTS_ASSERT(exception_thrown);
    };
    check5(0, "invalid number of components");
#endif /* MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION */
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
