/*!
 * \file   tests/AlgorithmsTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   22/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <array>
#include <memory>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Function/Algorithms.hxx"

struct AlgorithmsTest final : public tfel::tests::TestCase {
  AlgorithmsTest()
      : tfel::tests::TestCase("MGIS/Function", "AlgorithmsTests") {
  }  // end of AlgorithmsTest
  tfel::tests::TestResult execute() override {
    this->test1<1>();
    this->test1<2>();
    this->test1<3>();
    this->test1<4>();
    this->test1<5>();
    this->test1<6>();
    this->test1<7>();
    this->test1<8>();
    this->test1<9>();
    this->test1<10>();
    this->test1<11>();
    this->test1<12>();
    return this->result;
  }

 private:
  template <mgis::size_type N>
  void test1(){
    using namespace mgis;
    auto in = std::array<size_type, N> {};
    auto out = std::array<size_type, N>{};
    std::iota(in.begin(), in.end(), N);
    function::algorithm::copy<N>(in.begin(), in.end(), out.begin());
    TFEL_TESTS_CHECK(in == out);
  }
};

TFEL_TESTS_GENERATE_PROXY(AlgorithmsTest, "AlgorithmsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("AlgorithmsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
