/*!
 * \file   SharedSpaceTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/06/2025
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
#include "MGIS/Function/SharedSpace.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"

struct SharedSpaceTest final : public tfel::tests::TestCase {
  SharedSpaceTest()
      : tfel::tests::TestCase("MGIS/Function", "SharedSpaceTests") {
  }  // end of SharedSpaceTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis::function;
    TFEL_TESTS_STATIC_ASSERT(SpaceConcept<SharedSpace<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        ElementSpaceConcept<SharedSpace<BasicLinearSpace>>);
    TFEL_TESTS_STATIC_ASSERT(
        LinearElementSpaceConcept<SharedSpace<BasicLinearSpace>>);
  }
  void test2() {
    using namespace mgis::function;
    TFEL_TESTS_STATIC_ASSERT(
        SpaceConcept<SharedSpace<BasicLinearQuadratureSpace<3>>>);
    TFEL_TESTS_STATIC_ASSERT(
        ElementSpaceConcept<SharedSpace<BasicLinearQuadratureSpace<3>>>);
    TFEL_TESTS_STATIC_ASSERT(
        LinearElementSpaceConcept<SharedSpace<BasicLinearQuadratureSpace<3>>>);
    TFEL_TESTS_STATIC_ASSERT(
        QuadratureSpaceConcept<SharedSpace<BasicLinearQuadratureSpace<3>>>);
    TFEL_TESTS_STATIC_ASSERT(LinearQuadratureSpaceConcept<
                             SharedSpace<BasicLinearQuadratureSpace<3>>>);
  }
  void test3() {
    using namespace mgis::function;
    auto qspace =
        SharedSpace{std::make_shared<BasicLinearQuadratureSpace<3>>(10)};
    TFEL_TESTS_CHECK_EQUAL(getSpaceSize(qspace), 30);
    TFEL_TESTS_CHECK_EQUAL(getNumberOfCells(qspace), 10);
    TFEL_TESTS_CHECK_EQUAL(getNumberOfQuadraturePoints(qspace, 2), 3);
    TFEL_TESTS_CHECK_EQUAL(getQuadraturePointOffset(qspace, 2, 1), 7);
  }
};

TFEL_TESTS_GENERATE_PROXY(SharedSpaceTest, "SharedSpaceTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("SharedSpaceTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
