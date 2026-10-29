/*!
 * \file   RotateFunctionsTest.cxx
 * \brief
 * \author Thomas Heler
 * \date   08/06/2020
 */

#include <cmath>
#include <array>
#include <numbers>
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
    // rotation of 30 degrees around z: r stores column-major the matrix
    // R = {{c, s, 0}, {-s, c, 0}, {0, 0, 1}}, which is not symmetric
    constexpr auto cste = std::numbers::sqrt2_v<real>;
    const auto c = std::cos(std::numbers::pi_v<real> / 6);
    const auto s = std::sin(std::numbers::pi_v<real> / 6);
    const std::array<real, 9> r = {c, -s, 0, s, c, 0, 0, 0, 1};
    // gradient with a shear component (stored as sqrt(2) times exy)
    const std::array<real, 4> ge = {1, 2, 3, 4};
    std::array<real, 4> me;
    rotateGradients(me, b, ge, r);
    // reference: R E R^T
    const auto exy = ge[3] / cste;
    assert_equal(me[0], c * c * ge[0] + 2 * c * s * exy + s * s * ge[1]);
    assert_equal(me[1], s * s * ge[0] - 2 * c * s * exy + c * c * ge[1]);
    assert_equal(me[2], ge[2]);
    assert_equal(me[3],
                 cste * (c * s * (ge[1] - ge[0]) + (c * c - s * s) * exy));
    // the rotation of the thermodynamic forces is the inverse rotation
    std::array<real, 4> gs;
    rotateThermodynamicForces(gs, b, me, r);
    for (size_type i = 0; i != 4; ++i) {
      assert_equal(gs[i], ge[i]);
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
