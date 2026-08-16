/*!
 * \file   QuantityFunctionTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   16/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
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
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/TFEL/Quantity.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"

struct QuantityFunctionsTest final : public tfel::tests::TestCase {
  QuantityFunctionsTest()
      : tfel::tests::TestCase("MGIS/Function", "QuantityFunctionsTests") {
  }  // end of QuantityFunctionsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
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
    using namespace mgis;
    using namespace mgis::function;
    using time = ::tfel::math::qt<::tfel::math::unit::Time, real>;
    TFEL_TESTS_STATIC_ASSERT(
        (std::same_as<
            evaluator_result<QuantityModifier<FunctionView<BasicLinearSpace>,
                                              ::tfel::math::unit::Time>>,
            ::tfel::math::const_qt_ref<::tfel::math::unit::Time, real>>));
    TFEL_TESTS_STATIC_ASSERT(
        (std::is_invocable_v<
            mgis::function::internals::MultiplyByScalarOperator,
            ::tfel::math::const_qt_ref<::tfel::math::unit::Time, real>>));
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    constexpr auto values = []() -> std::optional<std::array<real, 4>> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{4};
      Function f(space, 1);
      Function f2(space, 1);
      auto t = f | as_qt<::tfel::math::unit::Time>;
      auto t2 = f2 | as_qt<::tfel::math::unit::Time>;
      t(0) = time{1};
      t(1) = time{-2};
      t(2) = time{-5};
      t(3) = time{4};
      const auto ok = assign(ctx, t2, t | multiply_by_scalar(2));
      if (!ok) {
        return {};
      }
      return std::array{f2(0)[0], f2(1)[0], f2(2)[0], f2(3)[0]};
    }();
    TFEL_TESTS_STATIC_ASSERT(values.has_value());
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[0], 2));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[1], -4));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[2], -10));
    TFEL_TESTS_STATIC_ASSERT(check_value((*values)[3], 8));
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    using stress = ::tfel::math::qt<::tfel::math::unit::Stress, real>;
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    constexpr auto values = []() -> std::optional<std::array<real, 4>> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{1};
      Function f(space, 4);
      Function f2(space, 4);
      auto t = f | as_stensor<2u, stress>;
      auto t2 = f2 | as_stensor<2u, stress>;
      t(0) = {stress{1}, stress{-2}, stress{-5}, stress{4}};
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
};

TFEL_TESTS_GENERATE_PROXY(QuantityFunctionsTest, "QuantityFunctionsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("QuantityFunctionsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}

