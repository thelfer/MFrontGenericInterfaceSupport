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
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

struct ImmutableFunctionTest final : public tfel::tests::TestCase {
  ImmutableFunctionTest()
      : tfel::tests::TestCase("MGIS/Function", "ImmutableFunctionTests") {
  }  // end of ImmutableFunctionTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    auto space = std::make_shared<BasicLinearSpace>(0);
    auto values = std::vector<real>{};
    auto f = ImmutableFunctionView<BasicLinearSpace, {}>(space, values);
  }
};

TFEL_TESTS_GENERATE_PROXY(ImmutableFunctionTest, "ImmutableFunctionTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FunctionTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
