/*!
 * \file   HDF5Test.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/01/2026
 */

#include <map>
#include <array>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <iostream>
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Utilities/HDF5Support.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  using namespace mgis::utilities::hdf5;
  if (argc != 2) {
    std::cerr << "HDF5Test: invalid number of arguments\n";
    std::exit(-1);
  }
  auto ctx = Context{};
  auto file = H5::H5File("HDF5Test.h5", H5F_ACC_TRUNC);
  auto root = file.openGroup("/");
  const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
  MaterialDataManager m{b, 100};
  const auto de = 5.e-5;
  // initialize the external state variable
  m.s1.external_state_variables["Temperature"] = 293.15;
  // copy d.s1 in d.s0
  update(m);
  const auto dt = real(180);
  for (size_type idx = 0; idx != m.n; ++idx) {
    m.s1.gradients[idx * m.s1.gradients_stride] = de;
  }
  for (size_type i = 0; i != 20; ++i) {
    auto og = createGroup(ctx, root, "TimeStep_" + std::to_string(i));
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
    for (size_type idx = 0; idx != m.n; ++idx) {
      m.s1.gradients[idx * m.s1.gradients_stride] += de;
    }
  }
  // restore and checks
  const auto o =
      getVariableOffset(b.isvs, "EquivalentViscoplasticStrain", b.hypothesis);
  auto pi = std::array<real, 20>{};  // values of the equivalent plastic strain
  auto pe = std::array<real, 20>{};  // values of the equivalent plastic strain
  const auto ni = size_type{o};
  const auto ne =
      size_type{(m.n - 1) * m.s0.internal_state_variables_stride + o};
  for (size_type i = 0; i != 20; ++i) {
    auto og = openGroup(ctx, root, "TimeStep_" + std::to_string(i));
    if (isInvalid(og)) {
      std::cerr << ctx.getErrorMessage() << '\n';
      return EXIT_FAILURE;
    }
    setExternalStateVariable(m.s1, "Temperature", 500);
    if (!restore(ctx, m.s1, *og, {.restore_mass_densities = false})) {
      std::cerr << ctx.getErrorMessage() << '\n';
      return EXIT_FAILURE;
    }
    const auto& T = m.s1.external_state_variables.at("Temperature");
    if (!std::holds_alternative<MaterialStateManager::MutableFieldHolder>(T)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto& Tvalue =
        std::get<MaterialStateManager::MutableFieldHolder>(T).value;
    if (!std::holds_alternative<real>(Tvalue)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    if (std::abs(std::get<real>(Tvalue) - 293.15) > 1e-10) {
      std::cerr << "invalid value for the temperature\n";
      return EXIT_FAILURE;
    }
    pi[i] = m.s1.internal_state_variables[ni];
    pe[i] = m.s1.internal_state_variables[ne];
  }
  //
  const auto p_ref = std::array<real, 20>{
      1.3523277308229e-11, 1.0955374667213e-07, 5.5890770166084e-06,
      3.2392193670428e-05, 6.645865307584e-05,  9.9676622883138e-05,
      0.00013302758358953, 0.00016635821069889, 0.00019969195920296,
      0.00023302522883648, 0.00026635857194317, 0.000299691903777,
      0.0003330252373404,  0.00036635857063843, 0.00039969190397718,
      0.00043302523730968, 0.00046635857064314, 0.00049969190397646,
      0.00053302523730979, 0.00056635857064313};
  std::cerr.precision(14);
  for (size_type i = 0; i != 20; ++i) {
    if (std::abs(pi[i] - p_ref[i]) > 1.e-12) {
      std::cerr << "invalid value for the equivalent "
                   "viscoplastic strain at the first integration point"
                << "(expected '" << p_ref[i] << "', computed '" << pi[i]
                << "')\n";
      return EXIT_FAILURE;
    }
    if (std::abs(pe[i] - p_ref[i]) > 1.e-12) {
      std::cerr << "invalid value for the equivalent "
                   "viscoplastic strain at the last integration point"
                << "(expected '" << p_ref[i] << "', computed '" << pe[i]
                << "')\n";
      return EXIT_FAILURE;
    }
  }
  // save and restore a MaterialDataManager
  auto m4 = MaterialDataManager{b, m.n};
  setExternalStateVariable(m4.s0, "Temperature", 0);
  setExternalStateVariable(m4.s1, "Temperature", 0);
  auto og4 = createGroup(ctx, root, "MaterialDataManager");
  if ((isInvalid(og4)) || (!save(ctx, *og4, m)) ||
      (!restore(ctx, m4, *og4, {{.restore_mass_densities = false}}))) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  if ((std::abs(m4.s0.internal_state_variables[ne] - p_ref[19]) > 1.e-12) ||
      (std::abs(m4.s1.internal_state_variables[ne] - p_ref[19]) > 1.e-12)) {
    std::cerr << "invalid values restored in a MaterialDataManager\n";
    return EXIT_FAILURE;
  }
  // external state variables which are not scalars, given at each
  // integration point
  const auto b2 = load(argv[1], "TensorialExternalStateVariableTest",
                       Hypothesis::TRIDIMENSIONAL);
  auto values = std::map<std::string, std::vector<real>>{};
  auto m2 = MaterialDataManager{b2, 2};
  auto m3 = MaterialDataManager{b2, 2};
  for (const auto& mp : b2.mps) {
    setMaterialProperty(m2.s1, mp.name, 1);
    setMaterialProperty(m3.s1, mp.name, 1);
  }
  for (const auto& v : b2.esvs) {
    auto& w = values[v.name];
    w.resize(m2.n * getVariableSize(v, b2.hypothesis));
    std::iota(w.begin(), w.end(), real{1});
    setExternalStateVariable(m2.s1, v.name, w);
    auto zeros = std::vector<real>(w.size(), real{0});
    setExternalStateVariable(m3.s1, v.name, zeros);
  }
  auto og = createGroup(ctx, root, "TensorialExternalStateVariables");
  if ((isInvalid(og)) || (!save(ctx, *og, m2.s1)) ||
      (!restore(ctx, m3.s1, *og, {.restore_mass_densities = false}))) {
    std::cerr << ctx.getErrorMessage() << '\n';
    return EXIT_FAILURE;
  }
  for (const auto& [n, w] : values) {
    const auto& esv = std::get<MaterialStateManager::MutableFieldHolder>(
        m3.s1.external_state_variables.at(n));
    if (std::get<std::vector<real>>(esv.value) != w) {
      std::cerr << "invalid values for the external state variable '" << n
                << "'\n";
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}  // end of main
