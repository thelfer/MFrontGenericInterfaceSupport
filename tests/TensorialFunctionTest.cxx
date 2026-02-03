/*!
 * \file   TensorialFunctionTest.cxx
 * \brief
 * \author th202608
 * \date   28/08/2025
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <array>
#include <cmath>
#include <memory>
#include <cstdlib>
#include <optional>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Material/Lame.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/FixedSizeView.hxx"
#include "MGIS/Function/FixedSizeModifier.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"
#include "MGIS/Function/TFEL/TensorialFunction.hxx"

struct TensorialFunctionsTest final : public tfel::tests::TestCase {
  TensorialFunctionsTest()
      : tfel::tests::TestCase("MGIS/Function", "TensorialFunctionsTests") {
  }  // end of TensorialFunctionsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    return this->result;
  }

 private:
  static constexpr bool check_value(const mgis::real a,
                                    const mgis::real b,
                                    const mgis::real eps = mgis::real{
                                        1e-12}) noexcept {
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    return local_abs(a - b) < eps;
  }  // end of check_value
  void test1() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto values = []() -> std::optional<std::array<real, 4>> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{1};
      Function f(space, 4);
      Function f2(space, 4);
      auto t = f | as_stensor<2>;
      auto t2 = f2 | as_stensor<2>;
      t(0) = {1, -2, -5, 4};
      const auto ok = assign(ctx, t2, t | multiply_by_scalar(2));
      if (!ok) {
        return {};
      }
      return std::array{f2(0)[0], f2(0)[1], f2(0)[2], f2(0)[3]};
    }();
    TFEL_TESTS_STATIC_ASSERT(values.has_value());
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[0], 2));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[1], -4));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[2], -10));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[3], 8));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test2() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    TFEL_TESTS_STATIC_ASSERT(
        (FunctionConcept<StensorFunction<BasicLinearSpace, 1u>>));
    constexpr auto values = []() {
      auto space = BasicLinearSpace{2};
      auto sig = StensorFunction<BasicLinearSpace, 1u>{space};
      sig(0) = tfel::math::stensor<1u, real>::Id();
      sig(1) = -tfel::math::stensor<1u, real>::Id();
      auto sig_values = std::array<real, 6>{};
      std::copy(sig.data().begin(), sig.data().end(), sig_values.begin());
      return sig_values;
    }();
    TFEL_TESTS_STATIC_ASSERT(check_value(values[0], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[1], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[2], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[3], -1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[4], -1));
    TFEL_TESTS_STATIC_ASSERT(check_value(values[5], -1));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test3() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto values = []() -> std::optional<std::array<real, 6>> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{2};
      auto sig = StensorFunction<BasicLinearSpace, 1u>{space};
      auto sig2 = StensorFunction<BasicLinearSpace, 1u>{space};
      sig(0) = tfel::math::stensor<1u, real>::Id();
      sig(1) = -tfel::math::stensor<1u, real>::Id();
      const auto ok = assign(ctx, sig2, view(sig));
      if (!ok) {
        return {};
      }
      auto sig2_values = std::array<real, 6>{};
      std::copy(sig2.data().begin(), sig2.data().end(), sig2_values.begin());
      return sig2_values;
    }();
    TFEL_TESTS_STATIC_ASSERT(values.has_value());
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[0], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[1], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[2], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[3], -1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[4], -1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[5], -1));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test4() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto values = []() -> std::optional<std::array<real, 3>> {
      auto space = BasicLinearSpace{1};
      auto sig = StensorFunction<BasicLinearSpace, 1u>{space};
      auto sig2_v = view(sig);
      sig2_v(0) = tfel::math::stensor<1u, real>::Id();
      auto sig_values = std::array<real, 3>{};
      std::copy(sig.data().begin(), sig.data().end(), sig_values.begin());
      return sig_values;
    }();
    TFEL_TESTS_STATIC_ASSERT(values.has_value());
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[0], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[1], 1));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[2], 1));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test5() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace tfel::math;
    using namespace tfel::material;
    using namespace mgis;
    using namespace mgis::function;
    using mgis::real;
    constexpr auto E = real{150e9};
    constexpr auto nu = 0.3;
    constexpr auto l = computeLambda(E, nu);
    constexpr auto m = computeMu(E, nu);
    constexpr auto k = l + (2 * m) / 3;
    constexpr auto eps = E * real{1e-14};
    constexpr auto values = [&]() -> std::optional<std::array<real, 16>> {
      using Stensor4 = st2tost2<2u, real>;
      auto space = BasicLinearSpace{1};
      auto D = ST2toST2Function<BasicLinearSpace, 2u>{space};
      D(0) = k * Stensor4::IxI() + 2 * m * Stensor4::K();
      auto D_values = std::array<real, 16>{};
      std::copy(D.data().begin(), D.data().end(), D_values.begin());
      return D_values;
    }();
    TFEL_TESTS_STATIC_ASSERT(values.has_value());
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[0], l + 2 * m, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[1], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[2], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[3], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[4], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[5], l + 2 * m, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[6], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[7], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[8], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[9], l, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[10], l + 2 * m, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[11], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[12], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[13], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[14], 0, eps));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[15], 2 * m, eps));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
};

TFEL_TESTS_GENERATE_PROXY(TensorialFunctionsTest, "TensorialFunctionsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("TensorialFunctionsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
