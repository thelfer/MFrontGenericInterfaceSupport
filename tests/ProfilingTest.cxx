/*!
 * \file   tests/ProfilingTest.cxx
 * \brief  Unit tests for the MGIS profiling system via Context
 * \author Julien Rigal
 * \date   2026
 */

#include <cstdlib>
#include <iostream>
#include <thread>
#include <chrono>

#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"

#include "MGIS/Context.hxx"
#include "MGIS/Profiling.hxx"
#include "MGIS/ProfilingData.hxx"

struct ProfilingTest final : public tfel::tests::TestCase {
  ProfilingTest()
      : tfel::tests::TestCase("MGIS/Profiling", "ProfilingSystemTests") {
  }  // end of ProfilingTest

  tfel::tests::TestResult execute() override {
    this->testToggle();
    this->testTreeHierarchy();
    this->testAccumulation();
    this->testDisabledState();
    this->testLocalForcedState();
    return this->result;
  }

 private:
  void testToggle() {
    mgis::Context ctx;

    // By default
    ctx.enableProfiling(false);
    TFEL_TESTS_ASSERT(!ctx.isProfilingEnabled());

    ctx.enableProfiling(true);
    TFEL_TESTS_ASSERT(ctx.isProfilingEnabled());
  }

  void testTreeHierarchy() {
    mgis::Context ctx;
    ctx.enableProfiling(true);

    // Simulate a nested structure
    {
      CatchTimeSection(ctx, "Parent");
      { CatchTimeSection(ctx, "Child1"); }
      { CatchTimeSection(ctx, "Child2"); }
    }

    const auto& root = ctx.getProfilingResultTree();

    // The root must have exactly 1 child ("Parent")
    TFEL_TESTS_ASSERT(root.children.size() == 1);

    const auto& parent = root.children[0];
    TFEL_TESTS_CHECK_EQUAL(parent->name, "Parent");
    TFEL_TESTS_CHECK_EQUAL(parent->calls, 1u);
    TFEL_TESTS_ASSERT(parent->time_in_seconds >= 0.0);

    // "Parent" must have exactly 2 children ("Child1" and "Child2")
    TFEL_TESTS_ASSERT(parent->children.size() == 2);
    TFEL_TESTS_CHECK_EQUAL(parent->children[0]->name, "Child1");
    TFEL_TESTS_CHECK_EQUAL(parent->children[1]->name, "Child2");
  }

  void testAccumulation() {
    mgis::Context ctx;
    ctx.enableProfiling(true);

    // Simulate a loop to check that "calls" is incremented
    {
      CatchTimeSection(ctx, "Loop");
      for (int i = 0; i < 5; ++i) {
        CatchTimeSection(ctx, "Body");
      }
    }

    const auto& root = ctx.getProfilingResultTree();
    TFEL_TESTS_ASSERT(root.children.size() == 1);

    const auto& loop = root.children[0];
    TFEL_TESTS_CHECK_EQUAL(loop->name, "Loop");
    TFEL_TESTS_CHECK_EQUAL(loop->calls, 1u);

    TFEL_TESTS_ASSERT(loop->children.size() == 1);
    const auto& body = loop->children[0];
    TFEL_TESTS_CHECK_EQUAL(body->name, "Body");
    TFEL_TESTS_CHECK_EQUAL(
        body->calls, 5u);  // If the loop ran five times, we have a problem
  }

  void testDisabledState() {
    mgis::Context ctx;
    ctx.enableProfiling(false);  // Profiling disabled

    { CatchTimeSection(ctx, "ShouldNotAppear"); }

    // The tree must be empty because the block was inactive
    const auto& root = ctx.getProfilingResultTree();
    TFEL_TESTS_ASSERT(root.children.empty());
  }

  void testLocalForcedState() {
    mgis::Context ctx;
    ctx.enableProfiling(false);  // Global profiling disabled

    {
      // Force local activation
      CatchLocalTimeSection(ctx, "ForcedSection", true);
    }

    // The node then must exist
    const auto& root = ctx.getProfilingResultTree();
    TFEL_TESTS_ASSERT(root.children.size() == 1);
    TFEL_TESTS_CHECK_EQUAL(root.children[0]->name, "ForcedSection");
  }
};

TFEL_TESTS_GENERATE_PROXY(ProfilingTest, "ProfilingTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ProfilingSystemTests.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}