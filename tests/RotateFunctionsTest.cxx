/*!
 * \file   RotateFunctionsTest.cxx
 * \brief    
 * \author Thomas Heler
 * \date   08/06/2020
 */

#include <cmath>
#include <array>
#include <stdexcept>
#include <iostream>
#include "MGIS/Behaviour/Behaviour.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  bool success = true;
  auto assert_equal = [&success](const real& a, const real b) {
    constexpr const auto e = real(1.e-14);
    if (std::abs(a - b) > e) {
      success = false;
    }
  };
  if (argc != 2) {
    std::cerr << "RotateFunctionsTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto b = load(argv[1], "OrthotropicElasticity",
                        Hypothesis::GENERALISEDPLANESTRAIN);
    const std::array<real, 9> r = {0, 1, 0, 1, 0, 0, 0, 0, 1};
    const std::array<real, 4> ge = {1, 0, 0, 0};
    std::array<real, 4> me;
    rotateGradients(me, b, ge, r);
    assert_equal(me[0], 0);
    assert_equal(me[1], 1);
    assert_equal(me[2], 0);
    assert_equal(me[3], 0);
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
