/*!
 * \file   tests/ThreadPoolTest.cxx
 * \brief  This file contains unit tests of the `ThreadPool` class
 * \date   27/10/2026
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/ThreadPool.hxx"

struct ThreadPoolTest final : public tfel::tests::TestCase {
  ThreadPoolTest()
      : tfel::tests::TestCase("MGIS", "ThreadPoolTests") {
  }  // end of ThreadPoolTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    TFEL_TESTS_CHECK_THROW(mgis::ThreadPool{0}, std::runtime_error);
  }  // end of test1
};

TFEL_TESTS_GENERATE_PROXY(ThreadPoolTest, "ThreadPoolTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ThreadPoolTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
