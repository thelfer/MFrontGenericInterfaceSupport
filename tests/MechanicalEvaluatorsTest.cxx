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
#include "MGIS/Function/FixedSizeView.hxx"
#include "MGIS/Function/FixedSizeModifier.hxx"
#include "MGIS/Function/Mechanics.hxx"

struct MechanicalEvaluatorsTest final : public tfel::tests::TestCase {
  MechanicalEvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function", "MechanicalEvaluatorsTests") {
  }  // end of MechanicalEvaluatorsTest
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
    auto space = BasicLinearSpace(1);
    auto values = std::vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(space, values, 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 3);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    const auto s = view<3>(f) | as_stensor<1>;
    const auto e = vmis(s);
    TFEL_TESTS_ASSERT(std::abs(e(0) - std::sqrt(3)) < eps);
    const auto e2 = hydrostatic_stress(s);
    TFEL_TESTS_ASSERT(std::abs(e2(0) - 2) < eps);
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto space = BasicLinearSpace(1);
    const auto values = std::vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(space, values, 3);
    auto seq = Function<BasicLinearSpace>(space, 1);
    TFEL_TESTS_ASSERT(seq.isScalar());
    const auto ok = f | as_stensor<1> | vmis | seq;
    TFEL_TESTS_ASSERT(ok);
    TFEL_TESTS_ASSERT(std::abs(seq(0)[0] - std::sqrt(3)) < eps);
    auto seq2 = Function<BasicLinearSpace>(space, 1);
    const auto ok2 = vmis(as_stensor<1>(view<3>(f))) | seq2;
    TFEL_TESTS_ASSERT(ok2);
    TFEL_TESTS_ASSERT(std::abs(seq2(0)[0] - std::sqrt(3)) < eps);
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    Context ctx;
    auto space = BasicLinearSpace(1);
    Function<BasicLinearSpace> pk1_function(space, 9);
    Function<BasicLinearSpace> F_function(space, 9);
    Function<BasicLinearSpace> sig_function(space, 6);
    auto pk1 = pk1_function | as_tensor<3>;
    auto F = F_function | as_tensor<3>;
    auto sig = sig_function | as_stensor<3>;
    TFEL_TESTS_STATIC_ASSERT(FunctionConcept<std::decay_t<decltype(pk1)>>);
    TFEL_TESTS_STATIC_ASSERT(EvaluatorConcept<std::decay_t<decltype(pk1)>>);
    TFEL_TESTS_STATIC_ASSERT(
        number_of_components<std::decay_t<decltype(pk1)>> == 9);
    TFEL_TESTS_ASSERT(pk1.check(ctx));
    TFEL_TESTS_CHECK_EQUAL(pk1.getNumberOfComponents(), 9);
    TFEL_TESTS_ASSERT(F.check(ctx));
    TFEL_TESTS_CHECK_EQUAL(F.getNumberOfComponents(), 9);
    pk1(0) = tfel::math::tensor<3u, real>::zero();
    F(0) = tfel::math::tensor<3u, real>::Id();
    // Cauchy stress in material frame
    const auto cauchy = pk1 | from_pk1_to_cauchy(F);
    TFEL_TESTS_STATIC_ASSERT(EvaluatorConcept<std::decay_t<decltype(cauchy)>>);
    TFEL_TESTS_STATIC_ASSERT(
        number_of_components<std::decay_t<decltype(cauchy)>> == 6);
    TFEL_TESTS_ASSERT(cauchy.check(ctx));
    TFEL_TESTS_CHECK_EQUAL(cauchy.getNumberOfComponents(), 6);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[0] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[1] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[2] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[3] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[4] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(cauchy(0)[5] - 0) < eps);
    // evaluator of the Cauchy stress in global frame
    const auto R = tfel::math::tmatrix<3, 3>{{0, 1, 0},  //
                                             {0, 0, 1},
                                             {1, 0, 0}};
    const auto cauchy_r = pk1 | from_pk1_to_cauchy(F) | rotate(R);
    TFEL_TESTS_STATIC_ASSERT(
        EvaluatorConcept<std::decay_t<decltype(cauchy_r)>>);
    TFEL_TESTS_ASSERT(cauchy_r.check(ctx));
    // evaluation of the Cauchy stress in global frame
    const auto ok = pk1 | from_pk1_to_cauchy(F) | rotate(R) | sig;
    TFEL_TESTS_ASSERT(ok);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[0] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[1] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[2] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[3] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[4] - 0) < eps);
    TFEL_TESTS_ASSERT(std::abs(sig(0)[5] - 0) < eps);
  }
  void test4() {
    using namespace mgis;
    using namespace mgis::function;
    Context ctx;
    auto space = BasicLinearSpace(1);
    const auto pk1_values = std::vector<real>{0, 0, 0, 0, 0, 0, 0, 0, 0};
    const auto F1_values = std::vector<real>{1, 1, 1, 0, 0, 0, 0, 0, 0};
    auto sig_values = std::vector<real>{0, 0, 0, 0, 0, 0, 0, 0, 0};
    const FunctionView<BasicLinearSpace, {}, false> pk1_function(space,
                                                                 pk1_values, 9);
    const FunctionView<BasicLinearSpace, {}, false> F_function(space, F1_values,
                                                               9);
    FunctionView<BasicLinearSpace> sig_function(space, sig_values, 6);
    const auto R = tfel::math::tmatrix<3, 3>{{0, 1, 0},  //
                                             {0, 0, 1},
                                             {1, 0, 0}};
    auto sig = sig_function | as_stensor<3>;
    const auto ok = pk1_function | as_tensor<3> |  //
                    from_pk1_to_cauchy(F_function | as_tensor<3>) |
                    rotate_backwards(R) | sig;
    TFEL_TESTS_ASSERT(ok);
  }
  void test5() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto size = []() constexpr {
      auto space = BasicLinearSpace{2};
      auto f = Function<BasicLinearSpace>{space, 3};
      auto v = f | as_tvector<3>;
      return v.getNumberOfComponents();
    }
    ();
    TFEL_TESTS_STATIC_ASSERT(size == 3);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
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
