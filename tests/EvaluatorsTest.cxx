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
    this->test4();
    this->test5();
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
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = std::make_shared<BasicLinearSpace>(2);
    auto values = std::vector<real>{1, 2, 3, 4};
    const auto f = ImmutableFunctionView<BasicLinearSpace>(space, values, 2) |
                    as_tvector<2>;
    TFEL_TESTS_ASSERT(f.check(ctx));
    auto values2 = std::vector<real>{1, -2, 3, -4};
    const auto f2 = ImmutableFunctionView<BasicLinearSpace>(space, values2, 2) |
                    as_tvector<2>;
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
  void test5() {
    //     using namespace mgis;
    //     using namespace mgis::function;
    //     constexpr auto eps = real{1e-14};
    //     Context ctx;
    //     auto space = std::make_shared<BasicLinearSpace>(1);
    //     auto values = std::vector<real>{1e-3, 2e-3, -5e-3, 4e-3};
    //     const auto strain =
    //         ImmutableFunctionView<BasicLinearSpace>(space, values, 4) |
    //         as_stensor<2>;
    //     auto values2 = std::vector<real>{0, 0, 0, 0};
    //     auto stress = FunctionView<BasicLinearSpace>(space, values2, 4) |
    //     as_stensor<2>; auto hooke_law = [l = 150e9,
    //                       m = 75e9](const tfel::math::StensorConcept auto& e)
    //                       {
    //       using namespace tfel::math;
    //       constexpr auto N = getSpaceDimension<decltype(e)>();
    //       constexpr auto id = stensor<N, real>{};
    //       return l * tfel::math::trace(e) * id + 2 * m * e;
    //     };
    //     assign(stress, strain | transform(hooke_law));
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
