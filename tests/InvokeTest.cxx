/*!
 * \file   tests/InvokeTest.cxx
 * \brief  This file contains unit test of the `invoke` function
 * \author Thomas Helfer
 * \date   20/02/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <cerrno>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Utilities/Invoke.hxx"

struct InvokeTest final : public tfel::tests::TestCase {
  InvokeTest()
      : tfel::tests::TestCase("MGIS/Utilities", "InvokeTests") {
  }  // end of InvokeTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    auto f = [](const int x) {
      if (x < 0) {
        raise<std::range_error>("argument shall be positive");
      }
      return 2 * x;
    };
    Context ctx;
    auto r1 = invoke(ctx, f, 12);
    TFEL_TESTS_ASSERT(r1.has_value());
    TFEL_TESTS_CHECK_EQUAL(*r1, 24);
    auto r2 = invoke(ctx, f, -12);
    TFEL_TESTS_ASSERT(!r2.has_value());
  }  // end of InvokeTest

  void test2() {
    using namespace mgis;
    auto f = [](const int x) {
      if (x < 0) {
        raise<std::range_error>("argument shall be positive");
      }
    };
    Context ctx;
    auto r1 = invoke(ctx, f, 12);
    TFEL_TESTS_ASSERT(r1);
    auto r2 = invoke(ctx, f, -12);
    TFEL_TESTS_ASSERT(!r2);
  }  // end of InvokeTest

  void test3() {
    using namespace mgis;
    Context ctx;
    const auto r1 =
        invokeCheckErrno(ctx, static_cast<double (*)(double)>(std::log), 1);
    TFEL_TESTS_ASSERT(r1.has_value());
    TFEL_TESTS_ASSERT(std::abs(*r1) < 1e-14);
    const auto r2 =
        invokeCheckErrno(ctx, static_cast<double (*)(double)>(std::log), -1);
    TFEL_TESTS_ASSERT(!r2.has_value());
    const auto r3 = MGIS_INVOKE_CHECK_ERRNO(
        ctx, static_cast<double (*)(double)>(std::log), 1);
    TFEL_TESTS_ASSERT(r3.has_value());
    TFEL_TESTS_ASSERT(std::abs(*r3) < 1e-14);
    const auto r4 = MGIS_INVOKE_CHECK_ERRNO(
        ctx, static_cast<double (*)(double)>(std::log), -1);
    TFEL_TESTS_ASSERT(!r4.has_value());
  }  // end of test3
};   // end of InvokeTest

TFEL_TESTS_GENERATE_PROXY(InvokeTest, "InvokeTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("InvokeTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
