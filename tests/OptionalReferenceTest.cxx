/*!
 * \file   tests/OptionalReferenceTest.cxx
 * \brief  This file contains unit test of the `OptionalReference` class
 * \author Thomas Helfer
 * \date   20/02/2026
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
#include "MGIS/Utilities/OptionalReference.hxx"

struct OptionalReferenceTest final : public tfel::tests::TestCase {
  OptionalReferenceTest()
      : tfel::tests::TestCase("MGIS/Utilities", "OptionalReferenceTests") {
  }  // end of OptionalReferenceTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    auto p = OptionalReference<int>{};
    TFEL_TESTS_CHECK_EQUAL(p, nullptr);
    TFEL_TESTS_CHECK_EQUAL(p.get(), nullptr);
    int a;
    p = OptionalReference<int>{&a};
    TFEL_TESTS_CHECK_EQUAL(p.get(), &a);
  }  // end of test1
};

TFEL_TESTS_GENERATE_PROXY(OptionalReferenceTest, "OptionalReferenceTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("OptionalReferenceTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
