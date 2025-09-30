/*!
 * \file   EvaluatorsTest.cxx
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
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/FixedSizeView.hxx"
#include "MGIS/Function/FixedSizeModifier.hxx"
#include "MGIS/Function/Tensors.hxx"

struct EvaluatorsTest final : public tfel::tests::TestCase {
  EvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function", "EvaluatorsTests") {
  }  // end of EvaluatorsTest
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
    this->test12();
    this->test13();
    this->test14();
    this->test15();
    this->test16();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{3};
    auto values = std::vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(space, values, 1);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 1);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 1);
    const auto e = FixedSizeModifier<FunctionEvaluator<BasicLinearSpace>, 1>(f);
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
        EvaluatorConcept<FunctionEvaluator<BasicLinearSpace>>);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{1};
    auto values = std::vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(space, values, 3);
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
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{2};
    auto values = std::vector<real>{1, 2, 3, 4};
    const auto f =
        FunctionEvaluator<BasicLinearSpace>(space, values, 2) | as_tvector<2>;
    TFEL_TESTS_ASSERT(f.check(ctx));
    auto values2 = std::vector<real>{1, -2, 3, -4};
    const auto f2 =
        FunctionEvaluator<BasicLinearSpace>(space, values2, 2) | as_tvector<2>;
    TFEL_TESTS_ASSERT(f2.check(ctx));
    //
    const auto n1 = f | negate;
    TFEL_TESTS_ASSERT(n1.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(n1(0)[0] + 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(n1(0)[1] + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(n1(1)[0] + 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(n1(1)[1] + 4) < eps);
    const auto s2 = f | multiply_by_scalar(4);
    TFEL_TESTS_ASSERT(s2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(s2(0)[0] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(s2(0)[1] - 8) < eps);
    TFEL_TESTS_ASSERT(std::abs(s2(1)[0] - 12) < eps);
    TFEL_TESTS_ASSERT(std::abs(s2(1)[1] - 16) < eps);
    const auto s3 = f | divide_by_scalar(3);
    TFEL_TESTS_ASSERT(s3.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(3 * s3(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(3 * s3(0)[1] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(3 * s3(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(3 * s3(1)[1] - 4) < eps);
    //
    const auto a1 = add(f, f2);
    TFEL_TESTS_ASSERT(a1.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(a1(0)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(a1(0)[1] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(a1(1)[0] - 6) < eps);
    TFEL_TESTS_ASSERT(std::abs(a1(1)[1] - 0) < eps);
    const auto a2 = f | add(f2);
    TFEL_TESTS_ASSERT(a2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(a2(0)[0] - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(a2(0)[1] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(a2(1)[0] - 6) < eps);
    TFEL_TESTS_ASSERT(std::abs(a2(1)[1] - 0) < eps);
    const auto s1 = f | substract(f2);
    TFEL_TESTS_ASSERT(s1.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(s1(0)[0] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(s1(0)[1] - 4) < eps);
    TFEL_TESTS_ASSERT(std::abs(s1(1)[0] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(s1(1)[1] - 8) < eps);
    const auto mv1 = f | mean_value(f2);
    TFEL_TESTS_ASSERT(mv1.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(mv1(0)[0] - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(mv1(0)[1] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(mv1(1)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(mv1(1)[1] - 0) < eps);
    const auto sp1 = f | inner_product(f2);
    TFEL_TESTS_ASSERT(sp1.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(sp1(0) + 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(sp1(1) + 7) < eps);
    const auto sp2 = inner_product(f, f2);
    TFEL_TESTS_ASSERT(sp2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(sp2(0) + 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(sp2(1) + 7) < eps);
    const auto n2 = f | divide(f | inner_product(f2));
    TFEL_TESTS_ASSERT(n2.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(3 * n2(0)[0] + 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(3 * n2(0)[1] + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(7 * n2(1)[0] + 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(7 * n2(1)[1] + 4) < eps);
    //
    const auto a3 = binary_operation(
        [](const tfel::math::VectorConcept auto& v1,
           const tfel::math::VectorConcept auto& v2) {
          return eval(v1 + 2 * v2);
        },
        f, f2);
    TFEL_TESTS_ASSERT(a3.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(a3(0)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(a3(0)[1] + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(a3(1)[0] - 9) < eps);
    TFEL_TESTS_ASSERT(std::abs(a3(1)[1] + 4) < eps);
    //
    const auto a4 = f | binary_operation(
                            [](const tfel::math::VectorConcept auto& v1,
                               const tfel::math::VectorConcept auto& v2) {
                              return eval(v1 + 2 * v2);
                            },
                            f2);
    TFEL_TESTS_ASSERT(a4.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(a4(0)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(a4(0)[1] + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(a4(1)[0] - 9) < eps);
    TFEL_TESTS_ASSERT(std::abs(a4(1)[1] + 4) < eps);
    //
    const auto a5 = binary_operation(
        [c = 2](const tfel::math::VectorConcept auto& v1,
                const tfel::math::VectorConcept auto& v2) {
          return eval(v1 + c * v2);
        },
        f, f2);
    TFEL_TESTS_ASSERT(a5.check(ctx));
    TFEL_TESTS_ASSERT(std::abs(a5(0)[0] - 3) < eps);
    TFEL_TESTS_ASSERT(std::abs(a5(0)[1] + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(a5(1)[0] - 9) < eps);
    TFEL_TESTS_ASSERT(std::abs(a5(1)[1] + 4) < eps);
  }
  template <typename T>
  static constexpr auto test5_shall_compile_test1 = requires(T& f3) {
    f3 | mgis::function::as_stensor<2>;
  };
  template <typename T>
  static constexpr auto test5_shall_compile_test2 = requires(const T& f3) {
    f3 | mgis::function::as_stensor<2>;
  };
  template <typename T>
  static constexpr auto test5_shall_not_compile_test1 = !requires(T & f3) {
    std::move(f3) | mgis::function::as_stensor<2>;
  };
  void test5() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{1};
    auto values = std::vector<real>{1e-3, 2e-3, -5e-3, 4e-3};
    static_assert(FunctionConcept<FunctionView<BasicLinearSpace>>);
    FunctionView<BasicLinearSpace> f(space, values, 4);
    auto strain = f | as_stensor<2>;
    TFEL_TESTS_ASSERT(strain.check(ctx));
    static_assert(std::same_as<decltype(f(0)), std::span<real>>);

    static_assert(std::same_as<decltype(f(0).data()), real*>);
    auto value = f.data(unsafe, 0);
    static_assert(std::same_as<decltype(value), real*>);
    strain(0) = 1e-3 * tfel::math::stensor<2, real>::Id();
    TFEL_TESTS_ASSERT(std::abs(values[0] - 1e-3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[1] - 1e-3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[2] - 1e-3) < eps);
    TFEL_TESTS_ASSERT(std::abs(values[3] - 0) < eps);
    //
    constexpr auto K = real{150e9};
    const auto stress = strain | multiply_by_scalar(K);
    TFEL_TESTS_ASSERT(std::abs(stress(0)[0] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(stress(0)[1] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(stress(0)[2] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(stress(0)[3] - 0) < K * eps);
    //
    auto values2 = std::vector<real>{1e-3, 2e-3, -5e-3, 4e-3};
    FunctionView<BasicLinearSpace> f2(space, values2, 4);
    auto stress_view = f2 | as_stensor<2>;
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<decltype(stress_view),
                      TensorView<FunctionView<BasicLinearSpace>,
                                 tfel::math::stensor<2, double>>>));
    const auto ok = assign(ctx, stress_view, strain | multiply_by_scalar(K));
    TFEL_TESTS_ASSERT(ok);
    TFEL_TESTS_ASSERT(std::abs(values2[0] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[1] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[2] - 150e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[3] - 0) < K * eps);
    const auto ok2 = strain | multiply_by_scalar(2 * K) | stress_view;
    TFEL_TESTS_ASSERT(ok2);
    TFEL_TESTS_ASSERT(std::abs(values2[0] - 300e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[1] - 300e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[2] - 300e6) < K * eps);
    TFEL_TESTS_ASSERT(std::abs(values2[3] - 0) < K * eps);
    //
    TFEL_TESTS_STATIC_ASSERT(FunctionConcept<FunctionView<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        !FunctionConcept<const FunctionView<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(FunctionConcept<Function<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        !FunctionConcept<const Function<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(EvaluatorConcept<FunctionView<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        test5_shall_compile_test1<FunctionView<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        test5_shall_compile_test2<FunctionView<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        test5_shall_not_compile_test1<Function<BasicLinearSpace>>);
  }
  void test6() {
    using namespace mgis;
    using namespace mgis::function;
#ifndef MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
    Context ctx;
    auto space = BasicLinearSpace{1};
    auto values = std::vector<real>{1e-3, 2e-3, -5e-3, 4e-3};
    static_assert(FunctionConcept<FunctionView<BasicLinearSpace>>);
    FunctionView<BasicLinearSpace> f(space, values, 4);
    TFEL_TESTS_ASSERT(
        !(TensorView<FunctionView<BasicLinearSpace>,
                     tfel::math::stensor<1>>::checkPreconditions(ctx, f)));
#endif /* MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING */
  }
  void test7() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{1};
    auto values = std::vector<real>{1, -2, -5, 4};
    FunctionEvaluator<BasicLinearSpace> f(space, values, 2);
    const auto max = view<2>(f) | maximum_component;
    TFEL_TESTS_ASSERT(max.check(ctx));
    TFEL_TESTS_ASSERT(max.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(max(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(max(1) - 4) < eps);
    const auto min = view<2>(f) | minimum_component;
    TFEL_TESTS_ASSERT(min.check(ctx));
    TFEL_TESTS_ASSERT(min.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(min(0) + 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(min(1) + 5) < eps);
    const auto abs_max = view<2>(f) | absolute_value | maximum_component;
    TFEL_TESTS_ASSERT(abs_max.check(ctx));
    TFEL_TESTS_ASSERT(abs_max.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(abs_max(0) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(abs_max(1) - 5) < eps);
    const auto abs_min = view<2>(f) | absolute_value | minimum_component;
    TFEL_TESTS_ASSERT(abs_min.check(ctx));
    TFEL_TESTS_ASSERT(abs_min.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(abs_min(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(abs_min(1) - 4) < eps);
  }
  void test8() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace{1};
    auto values = std::vector<real>{1, -2, -5, 4};
    FunctionEvaluator<BasicLinearSpace> f(space, values, 4);
    const auto max_vp = f | as_stensor<2> | eigen_values<> | maximum_component;
    TFEL_TESTS_ASSERT(max_vp.check(ctx));
    TFEL_TESTS_ASSERT(max_vp.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(max_vp(0) - (std::sqrt(41) - 1) / 2) < eps);
    const auto min_vp = f | as_stensor<2> | eigen_values<> | minimum_component;
    TFEL_TESTS_ASSERT(min_vp.check(ctx));
    TFEL_TESTS_ASSERT(min_vp.getNumberOfComponents() == 1);
    TFEL_TESTS_ASSERT(std::abs(min_vp(0) + 5) < eps);
  }
  void test9() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    const auto values = std::vector<real>{1, -2, -5, 4};
    {
      auto space = BasicLinearSpace{4};
      FunctionEvaluator<BasicLinearSpace> f(space, values, 1);
      const auto abs_values = view<1>(f) | absolute_value;
      TFEL_TESTS_ASSERT(abs_values.check(ctx));
      TFEL_TESTS_ASSERT(abs_values.getNumberOfComponents() == 1);
      TFEL_TESTS_ASSERT(std::abs(abs_values(0) - 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(1) - 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(2) - 5) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(3) - 4) < eps);
    }
    {
      auto space = BasicLinearSpace{2};
      FunctionEvaluator<BasicLinearSpace> f(space, values, 2);
      const auto abs_values = view<2>(f) | absolute_value;
      TFEL_TESTS_ASSERT(abs_values.check(ctx));
      TFEL_TESTS_ASSERT(getNumberOfComponents(abs_values) == 2);
      TFEL_TESTS_ASSERT(std::abs(abs_values(0)[0] - 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(0)[1] - 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(1)[0] - 5) < eps);
      TFEL_TESTS_ASSERT(std::abs(abs_values(1)[1] - 4) < eps);
    }
  }
  void test10() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    TFEL_TESTS_STATIC_ASSERT(
        (EvaluatorConcept<FixedSizeView<Function<BasicLinearSpace>, 1>>));
    Function<BasicLinearSpace> f(BasicLinearSpace{4}, 1);
    f(0)[0] = 1;
    f(1)[0] = -2;
    f(2)[0] = -5;
    f(3)[0] = 4;
    const auto f_1 = view<1>(f);
    const auto abs_values = f_1 | absolute_value;
    TFEL_TESTS_ASSERT(std::abs(abs_values(0) - 1) < eps);
    TFEL_TESTS_ASSERT(std::abs(abs_values(1) - 2) < eps);
    TFEL_TESTS_ASSERT(std::abs(abs_values(2) - 5) < eps);
    TFEL_TESTS_ASSERT(std::abs(abs_values(3) - 4) < eps);
  }
  static constexpr bool check_value(const mgis::real a,
                                    const mgis::real b) noexcept {
    constexpr auto eps = mgis::real{1e-12};
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    return local_abs(a - b) < eps;
  }  // end of check_value
  void test11() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto values = []() constexpr {
      Function f(BasicLinearSpace{4}, 1);
      Function f2(BasicLinearSpace{1}, 4);
      Function f3(BasicLinearSpace{1}, 4);
      tfel::math::map<tfel::math::tvector<4, real>>(&f(0)[0]) = {1, -2, -5, 4};
      (f2 | as_stensor<2>)(0) = {1, -2, -5, 4};
      (f3 | as_stensor<2>)(0) = {-4, 5, -3, -2};
      const auto abs_values = view<1>(f) | absolute_value;
      const auto trace_values = f2 | as_stensor<2> | trace;
      const auto trace_values2 = f2 | as_stensor<2> | trace | negate;
      const auto trace_values3 = f2 | as_stensor<2> | trace | negate |
                                 multiply_by_scalar(2) | divide_by_scalar(3);
      const auto min_values = view<4>(f2) | minimum_component;
      const auto max_values = view<4>(f2) | maximum_component;
      const auto transform_values =
          view<1>(f) | transform([](const real& x) { return x * x / 2; });
      const auto add_values = add(f2 | as_stensor<2>, f3 | as_stensor<2>);
      return std::array<real, 14>{abs_values(0),    abs_values(1),     //
                                  abs_values(2),    abs_values(3),     //
                                  trace_values(0),  trace_values2(0),  //
                                  trace_values3(0), min_values(0),     //
                                  max_values(0),    transform_values(3),
                                  add_values(0)[0], add_values(0)[1],
                                  add_values(0)[2], add_values(0)[3]};
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(values[0], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[1], 2));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[2], 5));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[3], 4));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[4], -6));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[5], 6));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[6], 4));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[7], -5));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[8], 4));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[9], 8));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[10], -3));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[11], 3));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[12], -8));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[13], 2));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test12() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto value = []() constexpr {
      Function f(BasicLinearSpace{1}, 4);
      (f | as_stensor<2>)(0) = {1, -2, -5, 4};
      const auto op = as_stensor<2> | trace | negate | multiply_by_scalar(2) |
                      divide_by_scalar(3);
      return (f.view() | op)(0);
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(value, 4));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test13() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto value = []() constexpr {
      Function f(BasicLinearSpace{2}, 16);
      auto D = f | as_st2tost2<2>;
      D(0) = tfel::math::st2tost2<2>::K();
      D(1) = tfel::math::st2tost2<2>::M();
      return (D(0) + D(1))(1, 2);
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(value, -real{5} / 6));
    constexpr auto n = []() constexpr {
      Function f(BasicLinearSpace{2}, 16);
      auto D = f | as_st2tost2<2>;
      return D.getNumberOfComponents();
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(n == 16);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test14() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto values = []() constexpr {
      auto eh = ContractViolationHandler{};
      auto space = BasicLinearSpace{2};
      Function<BasicLinearSpace, 1> f1(space);
      if (!f1.fill(eh, {12, 13})) {
        throw(std::runtime_error("fill failed"));
      }
      Function<BasicLinearSpace, 1> f2(space);
      if (!f2.fill(eh, {1, -3})) {
        throw(std::runtime_error("fill failed"));
      }
      const auto a = add(view(f1), view(f2));
      return std::array{a(0), a(1)};
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(check_value(values[0], 13));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[1], 10));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test15() {
#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS
    using namespace mgis;
    using namespace mgis::function;
    auto ctx = Context{};
    const auto space = BasicLinearSpace{200};
    Function<BasicLinearSpace, 1> f1(space);
    Function<BasicLinearSpace, 1> f2(space);
    for (size_type i = 0; i != space.size(); ++i) {
      f1(i) = 2 * i;
    }
    const auto e = view(f1) | multiply_by_scalar(0.5);
    const auto ok = assign(ctx, f2, std::execution::par, e);
    TFEL_TESTS_ASSERT(ok);
    for (size_type i = 0; i != space.size(); ++i) {
      TFEL_TESTS_ASSERT(std::abs(f2(i) - i) < 1e-12);
    }
#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */
  }
  void test16() {
#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS
    using namespace mgis;
    using namespace mgis::function;
    auto ctx = Context{};
    const auto space = BasicLinearSpace{200};
    Function<BasicLinearSpace, 1> f(space);
    for (size_type i = 0; i != space.size(); ++i) {
      const auto v = static_cast<real>(i);
      f(i) = 14 - 2 * (v - 12) * (v - 12);
    }
    const auto ov = scalar_reduce(
        ctx, std::execution::par, view(f),
        [](const real a, const real b) { return std::max(a, b); },
        std::numeric_limits<real>::min());
    TFEL_TESTS_ASSERT(isValid(ov));
    if (isValid(ov)) {
      TFEL_TESTS_ASSERT(std::abs(*ov - 14) < 1e-13);
    }
#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */
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
