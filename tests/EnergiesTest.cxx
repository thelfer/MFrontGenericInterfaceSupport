/*!
 * \file   tests/EnergiesTest.cxx
 * \brief  tests of the stored and dissipated energies computed by a
 * behaviour, through a behaviour data and a material data manager
 * \date   03/11/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cmath>
#include <vector>
#include <optional>
#include <cstdlib>
#include <iostream>
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

static bool success = true;

static void check(const bool b, const char* const msg) {
  if (!b) {
    std::cerr << "EnergiesTest: " << msg << '\n';
    success = false;
  }
}  // end of check

static bool are_equal(const mgis::real a, const mgis::real b) {
  return std::abs(a - b) < 1.e-12 * (1 + std::abs(b));
}  // end of are_equal

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "EnergiesTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto b = load(argv[1], "EnergiesTest", Hypothesis::TRIDIMENSIONAL);
    check(b.computesStoredEnergy, "the stored energy shall be computed");
    check(b.computesDissipatedEnergy,
          "the dissipated energy shall be computed");
    // behaviour data
    BehaviourData d(b);
    for (auto* const s : {&d.s0, &d.s1}) {
      setMaterialProperty(*s, "k", 2);
      setExternalStateVariable(*s, "Temperature", 293.15);
    }
    d.dt = 1;
    for (int step = 1; step != 3; ++step) {
      const auto e = 0.1 * step;
      d.s1.gradients[0] = e;
      d.K[0] = 0;
      auto v = make_view(d);
      check(integrate(v, b) == 1, "integration failed");
      check(are_equal(d.s1.stored_energy, e * e),
            "invalid stored energy (behaviour data)");
      check(are_equal(d.s1.dissipated_energy, 200 * step),
            "invalid dissipated energy (behaviour data)");
      update(d);
      check(are_equal(d.s0.stored_energy, e * e) &&
                are_equal(d.s0.dissipated_energy, 200 * step),
            "invalid energies after update (behaviour data)");
    }
    // material data managers, with the energies at the end of the time step
    // stored locally or externally
    const auto n = size_type{3};
    auto k = std::vector<real>{1, 2, 3};
    auto stored_energies = std::vector<real>(n, -1);
    auto dissipated_energies = std::vector<real>(n, -1);
    auto i1 = MaterialDataManagerInitializer{};
    i1.s1.stored_energies = stored_energies;
    i1.s1.dissipated_energies = dissipated_energies;
    ThreadPool p(2);
    for (const auto external : {false, true}) {
      // a material data manager can be neither copied nor moved
      auto om = std::optional<MaterialDataManager>{};
      if (external) {
        om.emplace(b, n, i1);
      } else {
        om.emplace(b, n);
      }
      auto& m = *om;
      for (auto* const s : {&m.s0, &m.s1}) {
        setMaterialProperty(*s, "k", k);
        setExternalStateVariable(*s, "Temperature", 293.15);
      }
      for (int step = 1; step != 3; ++step) {
        for (size_type i = 0; i != n; ++i) {
          m.s1.gradients[i * m.s1.gradients_stride] = 0.1 * step * (i + 1);
        }
        const auto it = IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
        const auto r =
            (step == 1) ? integrate(m, it, 1, 0, n) : integrate(p, m, it, 1);
        check(r == 1, "integration failed");
        for (size_type i = 0; i != n; ++i) {
          const auto e = 0.1 * step * (i + 1);
          check(are_equal(m.s1.stored_energies[i], k[i] * e * e / 2),
                "invalid stored energy");
          check(are_equal(m.s1.dissipated_energies[i], 100 * k[i] * step),
                "invalid dissipated energy");
          if (external) {
            check(are_equal(stored_energies[i], k[i] * e * e / 2) &&
                      are_equal(dissipated_energies[i], 100 * k[i] * step),
                  "energies not stored in the external arrays");
          }
        }
        update(m);
        for (size_type i = 0; i != n; ++i) {
          check(are_equal(m.s0.stored_energies[i], m.s1.stored_energies[i]) &&
                    are_equal(m.s0.dissipated_energies[i],
                              m.s1.dissipated_energies[i]),
                "invalid energies after update");
        }
      }
      for (size_type i = 0; i != n; ++i) {
        m.s1.stored_energies[i] = -1;
        m.s1.dissipated_energies[i] = -1;
      }
      revert(m);
      for (size_type i = 0; i != n; ++i) {
        check(are_equal(m.s1.stored_energies[i], m.s0.stored_energies[i]) &&
                  are_equal(m.s1.dissipated_energies[i],
                            m.s0.dissipated_energies[i]),
              "invalid energies after revert");
      }
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
