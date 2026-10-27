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
    this->test2();
    this->test3();
    return this->result;
  }

 private:
  void test1() {
    TFEL_TESTS_CHECK_THROW(mgis::ThreadPool{0}, std::runtime_error);
  }  // end of test1
  void test2() {
    struct Task {
      const Task* operator()() const { return this; }
    };
    mgis::ThreadPool p{1};
    const auto task = Task{};
    auto r = p.addTask(task).get();
    // the task is executed on a copy of the given object
    TFEL_TESTS_ASSERT(*r != &task);
  }  // end of test2
  void test3() {
    mgis::ThreadPool p{1};
    auto r = p.addTask([] { return 42.5; }).get();
    auto r2 = r;
    TFEL_TESTS_CHECK_EQUAL(*r2, 42.5);
    auto e =
        p.addTask([]() -> double { throw std::runtime_error("failed"); }).get();
    auto e2 = e;
    TFEL_TESTS_ASSERT(!e2);
  }  // end of test3
};

TFEL_TESTS_GENERATE_PROXY(ThreadPoolTest, "ThreadPoolTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ThreadPoolTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
