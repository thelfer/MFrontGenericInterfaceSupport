/*!
 * \file   HDF5Test.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/01/2026
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  using namespace mgis::utilities::hdf5;
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  auto ctx = Context{};
  auto f = H5::H5File("HDF5Test.h5", H5F_ACC_TRUNC);
  auto r = f.openGroup("/");
  const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
  MaterialDataManager m{b, 100};
  const auto de = 5.e-5;
  // initialize the external state variable
  m.s1.external_state_variables["Temperature"] = 293.15;
  // copy d.s1 in d.s0
  update(m);
  for (size_type idx = 0; idx != m.n; ++idx) {
    m.s1.gradients[idx * m.s1.gradients_stride] = de;
  }
  const auto dt = real(180);
  for (size_type i = 0; i != 20; ++i) {
    auto og = createGroup(ctx, r, "TimeStep_" + std::string{i});
    if (isInvalid(og)) {
      std::cerr << ctx.getErrorMessage() << '\n';
      return EXIT_FAILURE;
    }
    integrate(m, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, dt, 0, m.n);
    if (!save(ctx, *og, m.s1)) {
      std::cerr << ctx.getErrorMessage() << '\n';
      return EXIT_FAILURE;
    }
    update(m);
  }
  return EXIT_SUCCESS;
}  // end of main
