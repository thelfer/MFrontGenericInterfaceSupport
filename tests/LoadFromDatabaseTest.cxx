/*!
 * \file   tests/LoadFromDatabaseTest.cxx
 * \brief  tests of the loading of behaviours from the database
 * \date   03/11/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cstdlib>
#include <iostream>
#include "MGIS/Context.hxx"
#include "MGIS/Database.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  bool success = true;
  auto check = [&success](const bool b, const std::string& msg) {
    if (!b) {
      success = false;
      std::cerr << msg << '\n';
    }
    return b;
  };
  if (!check(argc == 2, "expected two arguments")) {
    return EXIT_FAILURE;
  }
  try {
    mgis::getDatabase().analyseLibrary(argv[1]);
    // small strain behaviour
    const auto b = loadFromDatabase(
        {.name = "Norton", .hypothesis = Hypothesis::PLANESTRAIN});
    check(b.behaviour == "Norton", "invalid behaviour name");
    check(b.hypothesis == Hypothesis::PLANESTRAIN, "invalid hypothesis");
    // finite strain behaviour
    auto o = FiniteStrainBehaviourOptions{};
    o.stress_measure = FiniteStrainBehaviourOptions::PK1;
    o.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
    const auto fb =
        loadFromDatabase(o, {.name = "FiniteStrainSingleCrystal",
                             .hypothesis = Hypothesis::TRIDIMENSIONAL});
    check(fb.behaviour == "FiniteStrainSingleCrystal",
          "invalid behaviour name");
    if (check(fb.thermodynamic_forces.size() == 1,
              "invalid number of thermodynamic forces")) {
      check(fb.thermodynamic_forces[0].name == "FirstPiolaKirchhoffStress",
            "invalid stress measure");
    }
#ifndef MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
    // a name matching several behaviours is rejected
    auto ctx = mgis::Context{};
    check(!loadFromDatabase(ctx, {.name = "Norton.*",
                                  .hypothesis = Hypothesis::TRIDIMENSIONAL})
               .has_value(),
          "a name matching several behaviours shall be rejected");
#endif /* MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING */
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
