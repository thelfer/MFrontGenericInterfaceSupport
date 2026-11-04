/*!
 * \file   tests/TangentOperatorTest.cxx
 * \brief  tests of the tangent and prediction operators returned by the
 * integrations of a material data manager and of a behaviour data, with the
 * `TangentKind` behaviour whose operator identifies the type of integration
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
#include <cstdlib>
#include <iostream>
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

static bool status = true;

static void check(const bool b, const char* const msg) {
  if (!b) {
    std::cerr << "TangentOperatorTest: " << msg << '\n';
    status = false;
  }
}  // end of check

static bool are_equal(const mgis::real a, const mgis::real b) {
  return std::abs(a - b) < 1.e-12 * (1 + std::abs(b));
}  // end of are_equal

/*!
 * \return the component (i, j) of the operator computed by the `TangentKind`
 * behaviour
 * \param[in] t: type of integration
 * \param[in] k: value of the material property `k`
 */
static mgis::real expected_operator(const mgis::behaviour::IntegrationType t,
                                    const mgis::real k,
                                    const mgis::size_type i,
                                    const mgis::size_type j) {
  const auto c = static_cast<int>(t);
  auto v = (i == j) ? c * k : mgis::real{0};
  if ((c > 0) && (i == 0) && (j == 1)) {
    v += k / 2;
  }
  return v;
}  // end of expected_operator

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "TangentOperatorTest: invalid number of arguments\n";
    std::exit(-1);
  }
  const IntegrationType types[] = {
      IntegrationType::PREDICTION_TANGENT_OPERATOR,
      IntegrationType::PREDICTION_SECANT_OPERATOR,
      IntegrationType::PREDICTION_ELASTIC_OPERATOR,
      IntegrationType::INTEGRATION_ELASTIC_OPERATOR,
      IntegrationType::INTEGRATION_SECANT_OPERATOR,
      IntegrationType::INTEGRATION_TANGENT_OPERATOR,
      IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR};
  try {
    const auto b = load(argv[1], "TangentKind", Hypothesis::TRIDIMENSIONAL);
    const auto ts = getStensorSize(b.hypothesis);
    // material data manager: the material property k differs at each point
    const auto n = size_type{5};
    auto k = std::vector<real>(n);
    for (size_type i = 0; i != n; ++i) {
      k[i] = 10 * (i + 1);
    }
    ThreadPool p(2);
    for (const auto t : types) {
      const auto prediction = static_cast<int>(t) < 0;
      for (const auto sos : {false, true}) {
        // full range, range [1, 4) and thread pool
        for (const auto mode : {0, 1, 2}) {
          MaterialDataManager m(b, n);
          for (auto* const s : {&m.s0, &m.s1}) {
            setMaterialProperty(*s, "k", k);
            setExternalStateVariable(*s, "Temperature", 293.15);
          }
          if (sos) {
            setMassDensity(m.s0, 2);
            setMassDensity(m.s1, 3);
          }
          for (size_type i = 0; i != n; ++i) {
            m.s1.gradients[i * m.s1.gradients_stride] = 1e-3 * (i + 1);
          }
          auto o = BehaviourIntegrationOptions{};
          o.integration_type = t;
          o.compute_speed_of_sound = sos;
          const auto ib = (mode == 1) ? size_type{1} : size_type{0};
          const auto ie = (mode == 1) ? size_type{4} : n;
          const auto r = (mode == 2) ? integrate(p, m, o, 1).exit_status
                                     : integrate(m, o, 1, ib, ie).exit_status;
          check(r == 1, "integration failed");
          check(m.K.size() == n * m.K_stride,
                "invalid size of the array of tangent operators");
          if (m.K.size() != n * m.K_stride) {
            continue;
          }
          for (size_type i = 0; i != n; ++i) {
            // the points outside the range are left untouched
            const auto treated = (i >= ib) && (i < ie);
            const auto K = m.K.data() + i * m.K_stride;
            for (size_type row = 0; row != ts; ++row) {
              for (size_type col = 0; col != ts; ++col) {
                const auto Kv =
                    treated ? expected_operator(t, k[i], row, col) : real{0};
                check(are_equal(K[row * ts + col], Kv),
                      "invalid value of the operator");
              }
            }
            // a prediction leaves the state at the end of the time step
            // untouched
            const auto integrated = treated && (!prediction);
            check(
                are_equal(m.s1.internal_state_variables[i], integrated ? 1 : 0),
                "invalid value of the internal state variable");
            check(are_equal(m.s1.thermodynamic_forces[i * ts],
                            integrated ? k[i] * 1e-3 * (i + 1) : 0),
                  "invalid value of the stress");
            if (sos) {
              // mass density at the beginning of the time step for a
              // prediction, at the end of the time step otherwise
              const auto v = prediction ? 2 * k[i] : 3 * k[i];
              check(are_equal(m.speed_of_sound[i], treated ? v : 0),
                    "invalid value of the speed of sound");
            }
          }
        }
      }
    }
    // behaviour data
    for (const auto t : types) {
      BehaviourData d(b);
      for (auto* const s : {&d.s0, &d.s1}) {
        setMaterialProperty(*s, "k", 7);
        setExternalStateVariable(*s, "Temperature", 293.15);
      }
      d.s1.gradients[0] = 1e-3;
      d.dt = 1;
      d.K[0] = static_cast<int>(t);
      auto v = make_view(d);
      check(integrate(v, b) == 1, "integration failed");
      for (size_type i = 0; i != ts; ++i) {
        for (size_type j = 0; j != ts; ++j) {
          check(are_equal(d.K[i * ts + j], expected_operator(t, 7, i, j)),
                "invalid value of the operator of the behaviour data");
        }
      }
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return status ? EXIT_SUCCESS : EXIT_FAILURE;
}
