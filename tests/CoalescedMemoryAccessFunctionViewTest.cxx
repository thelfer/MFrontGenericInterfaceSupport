/*!
 * \file   tests/CoalescedMemoryAccessFunctionViewTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   26/10/2025
 */

#include <cmath>
#include <memory>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Material/Lame.hxx"
#include "MGIS/Function/SharedSpace.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/Tensors.hxx"

namespace mgis::function {}  // end of namespace mgis::function

struct CoalescedMemoryAccessFunctionViewBaseTest final
    : public tfel::tests::TestCase {
  CoalescedMemoryAccessFunctionViewBaseTest()
      : tfel::tests::TestCase("MGIS/Function",
                              "CoalescedMemoryAccessFunctionViewBaseTests") {
  }  // end of CoalescedMemoryAccessFunctionViewBaseTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r =
        []() -> std::tuple<std::array<real, 2>, std::array<real, 2>> {
      auto space = BasicLinearSpace{2};
      auto c0 = Function<BasicLinearSpace, 1>{space};
      auto c1 = Function<BasicLinearSpace, 1>{space};
      c0(0) = 5;
      c0(1) = 12;
      c1(0) = -2;
      c1(1) = 3;
      auto components = std::array{view(c0), view(c1)};
      auto coalesced_view =
          CoalescedMemoryAccessFunctionViewBase<BasicLinearSpace, 2, true>(
              components);
      auto ptr0 = coalesced_view.getValuesPointers(0);
      auto ptr1 = coalesced_view.getValuesPointers(1);
      return {std::array{*ptr0[0], *ptr0[1]},  //
              std::array{*ptr1[0], *ptr1[1]}};
    }();
    TFEL_TESTS_STATIC_ASSERT(local_abs(std::get<0>(r)[0] - 5) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(std::get<0>(r)[1] + 2) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(std::get<1>(r)[0] - 12) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(std::get<1>(r)[1] - 3) < 1e-14);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
};

struct CoalescedMemoryAccessTensorViewTest final
    : public tfel::tests::TestCase {
  CoalescedMemoryAccessTensorViewTest()
      : tfel::tests::TestCase("MGIS/Function",
                              "CoalescedMemoryAccessTensorViewTests") {
  }  // end of CoalescedMemoryAccessTensorViewTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r = []() -> std::array<real, 3u> {
      auto ctx = ContractViolationHandler{};
      const auto ne = size_type{2};
      auto space = BasicLinearSpace{ne};
      std::array<const real, 8> values = {1, 10, 2, 20, 3, 30, 4, 40};
      const auto oscalar_functions =
          splitArrayIntoScalarFunctionViews<4>(ctx, space, values);
      if (!oscalar_functions.has_value()) {
        raise("splitArrayIntoScalarFunctionViews failed");
      }
      const auto f =
          CoalescedMemoryAccessTensorView<BasicLinearSpace,
                                          tfel::math::stensor<2, real>, false>{
              *oscalar_functions};
      return std::array{tfel::math::trace(f(0)), tfel::math::trace(f(1)),
                        tfel::math::trace(f(1) - 3 * f(0))};
    }();
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 6) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] - 60) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2] - 42) < 1e-14);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
};

struct CoalescedMemoryAccessCompositeTensorsViewTest final
    : public tfel::tests::TestCase {
  CoalescedMemoryAccessCompositeTensorsViewTest()
      : tfel::tests::TestCase(
            "MGIS/Function", "CoalescedMemoryAccessCompositeTensorsViewTests") {
  }  // end of CoalescedMemoryAccessCompositeTensorsViewTest
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
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r = []() -> std::array<real, 3u> {
      auto ctx = ContractViolationHandler{};
      const auto ne = size_type{2};
      auto space = BasicLinearSpace{ne};
      std::array<const real, 8> values = {1, 10, 2, 20, 3, 30, 4, 40};
      const auto oscalar_functions =
          splitArrayIntoScalarFunctionViews<4>(ctx, space, values);
      auto f =
          CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4, false>{
              *oscalar_functions};
      return std::array{
          tfel::math::trace(f.get<0, tfel::math::stensor<2, real>>(0)),
          tfel::math::trace(f.get<0, tfel::math::stensor<2, real>>(1)),
          tfel::math::trace(f.get<0, tfel::math::stensor<2, real>>(1) -
                            3 * f.get<0, tfel::math::stensor<2, real>>(0))};
    }();
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 6) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] - 60) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2] - 42) < 1e-14);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test2() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r = []() -> std::array<real, 8u> {
      auto ctx = ContractViolationHandler{};
      const auto ne = size_type{2};
      auto space = BasicLinearSpace{ne};
      std::array<real, 8> out_values = {0, 0, 0, 0, 0, 0, 0, 0};
      std::array<const real, 8> in_values = {1, 10, 2, 20, 3, 30, 4, 40};
      auto out = CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4>{
          space, out_values};
      const auto in =
          CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4, false>{
              space, in_values};
      for (size_type idx = 0; idx != ne; ++idx) {
        auto o = out.get<0, tfel::math::stensor<2, real>>(idx);
        const auto i = in.get<0, tfel::math::stensor<2, real>>(idx);
        o = i;
      }
      return out_values;
    }();
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 1) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] - 10) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2] - 2) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[3] - 20) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[4] - 3) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[5] - 30) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[6] - 4) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[7] - 40) < 1e-14);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
  void test3() {
    using namespace mgis;
    using namespace mgis::function;
    using CompositeView =
        CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4>;
    auto space = BasicLinearSpace{2};
    std::array<real, 6> out_values = {0, 0, 0, 0, 0, 0};
    auto ctx = Context{};
    TFEL_TESTS_ASSERT(
        !CompositeView::checkPreconditions(ctx, space, out_values));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "invalid number of values");
    if constexpr (config::contract_violation_policy ==
                  config::ContractViolationPolicy::RAISE) {
      bool has_thrown = false;
      try {
        auto out =
            CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 4>{
                space, out_values};
      } catch (std::exception&) {
        has_thrown = true;
      }
      TFEL_TESTS_ASSERT(has_thrown);
    }
  }  // end of test3
  void test4() {
    // elasticity
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    using CompositeView =
        CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 6>;
    using ImmutableCompositeView =
        CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 6, false>;
    constexpr auto young = real{150e9};
    constexpr auto nu = real{1} / 3;
    constexpr auto lambda = tfel::material::computeLambda(young, nu);
    constexpr auto mu = tfel::material::computeMu(young, nu);
    constexpr auto exx = 1e-3;
    constexpr auto eyy = -real{1e-3} / 3;
    constexpr auto seps = young * 1e-14;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r = []() -> std::array<real, 12u> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{2};
      auto sig_values = std::array<real, 12>{0, 0, 0, 0, 0, 0,  //
                                             0, 0, 0, 0, 0, 0};
      const auto eto_values = std::array<const real, 12>{exx, -exx,  //
                                                         eyy, -exx,  //
                                                         eyy, -exx,  //
                                                         0,   0,     //
                                                         0,   0,     //
                                                         0,   0};
      auto eto_view = ImmutableCompositeView{space, eto_values};
      auto sig_view = CompositeView{space, sig_values};
      for (size_type idx = 0; idx != getSpaceSize(space); ++idx) {
        using Stensor = tfel::math::stensor<3u, real>;
        constexpr auto id = Stensor::Id();
        const auto eto = eto_view.get<0, Stensor>(idx);
        auto sig = sig_view.get<0, Stensor>(idx);
        sig = lambda * tfel::math::trace(eto) * id + 2 * mu * eto;
      }
      return sig_values;
    }();
    constexpr auto K = (3 * lambda + 2 * mu) / 3;
    constexpr auto pr = K * (-3 * exx);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 150e6) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[4]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[6]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[8]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[10]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] - pr) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[3] - pr) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[5] - pr) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[7]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[9]) < seps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[11]) < seps);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }    // end of test4
  void test5() {
#ifndef MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS
    using namespace mgis;
    using namespace mgis::function;
    using ImmutableCompositeView =
        CoalescedMemoryAccessCompositeTensorsView<BasicLinearSpace, 2, false>;
    auto local_abs = [](const mgis::real r) { return r > 0 ? r : -r; };
    constexpr auto r = []() -> std::array<real, 4u> {
      auto ctx = ContractViolationHandler{};
      auto space = BasicLinearSpace{2};
      const auto values = std::array<const real, 4>{2, 3, -1, 4};
      const auto view = ImmutableCompositeView{space, values};
      return std::array{view.get<0, real>(0), view.get<1, real>(0),
                        view.get<0, real>(1), view.get<1, real>(1)};
    }();
    constexpr auto eps = real{1e-14};
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 2) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] + 1) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2] - 3) < eps);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[3] - 4) < eps);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }    // end of test5
};

TFEL_TESTS_GENERATE_PROXY(CoalescedMemoryAccessFunctionViewBaseTest,
                          "CoalescedMemoryAccessFunctionViewBaseTest");
TFEL_TESTS_GENERATE_PROXY(CoalescedMemoryAccessTensorViewTest,
                          "CoalescedMemoryAccessTensorViewTest");
TFEL_TESTS_GENERATE_PROXY(CoalescedMemoryAccessCompositeTensorsViewTest,
                          "CoalescedMemoryAccessCompositeTensorsViewTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("CoalescedMemoryAccessFunctionTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
