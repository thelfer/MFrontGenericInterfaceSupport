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
#include "MGIS/Function/SharedSpace.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"
#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"

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
          CoalescedMemoryAccessFunctionViewBase<BasicLinearSpace, 2>(
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
      using ScalarFunctionView =
          FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{
                                             .data_size = 1, .data_stride = 1}>;
      const auto ne = size_type{2};
      auto space = BasicLinearSpace{ne};
      real values[8] = {1, 10, 2, 20, 3, 30, 4, 40};
      auto c0 = ScalarFunctionView{space, std::span{values, ne}};
      auto c1 = ScalarFunctionView{space, std::span{values + ne, ne}};
      auto c2 = ScalarFunctionView{space, std::span{values + 2 * ne, ne}};
      auto c3 = ScalarFunctionView{space, std::span{values + 3 * ne, ne}};
      auto f = CoalescedMemoryAccessTensorView<BasicLinearSpace,
                                               tfel::math::stensor<2, real>>{
          std::array{c0, c1, c2, c3}};
      return std::array{tfel::math::trace(f(0)), tfel::math::trace(f(1)),
                        tfel::math::trace(f(1) - 3 * f(0))};
    }();
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[0] - 6) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[1] - 60) < 1e-14);
    TFEL_TESTS_STATIC_ASSERT(local_abs(r[2] - 42) < 1e-14);
#endif /* MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS */
  }
};

TFEL_TESTS_GENERATE_PROXY(CoalescedMemoryAccessFunctionViewBaseTest,
                          "CoalescedMemoryAccessFunctionViewBaseTest");
TFEL_TESTS_GENERATE_PROXY(CoalescedMemoryAccessTensorViewTest,
                          "CoalescedMemoryAccessTensorViewTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("CoalescedMemoryAccessFunctionTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
