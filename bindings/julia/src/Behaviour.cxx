/*!
 * \file   bindings/julia/src/Behaviour.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <jlcxx/jlcxx.hpp>
#include <jlcxx/const_array.hpp>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareBehaviour(jlcxx::Module &m) {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::FiniteStrainBehaviourOptions;
  using mgis::behaviour::Hypothesis;
  //
  Behaviour (*load)(const std::string &, const std::string &,
                    const Hypothesis) = &mgis::behaviour::load;
  Behaviour (*load2)(const FiniteStrainBehaviourOptions &, const std::string &,
                     const std::string &, const Hypothesis) =
      &mgis::behaviour::load;
  void (*setParameter1)(const Behaviour &, const std::string &, const double) =
      &mgis::behaviour::setParameter;
  void (*setParameter2)(const Behaviour &, const std::string &, const int) =
      &mgis::behaviour::setParameter;
  double (*getParameterDefaultValue1)(const Behaviour &, const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<double>;
  int (*getParameterDefaultValue2)(const Behaviour &, const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<int>;
  unsigned short (*getParameterDefaultValue3)(const Behaviour &,
                                              const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<unsigned short>;
  bool (*hasBounds)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasBounds;
  bool (*hasLowerBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasLowerBound;
  bool (*hasUpperBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasUpperBound;
  long double (*getLowerBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::getLowerBound;
  long double (*getUpperBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::getUpperBound;
  bool (*hasPhysicalBounds)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasPhysicalBounds;
  bool (*hasLowerPhysicalBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasLowerPhysicalBound;
  bool (*hasUpperPhysicalBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::hasUpperPhysicalBound;
  long double (*getLowerPhysicalBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::getLowerPhysicalBound;
  long double (*getUpperPhysicalBound)(const Behaviour &, const std::string &) =
      mgis::behaviour::getUpperPhysicalBound;
  //!
  m.add_bits<FiniteStrainBehaviourOptions::StressMeasure>(
      "FiniteStrainBehaviourOptionsStressMeasure");
  m.set_const("Cauchy", FiniteStrainBehaviourOptions::CAUCHY);
  m.set_const("CAUCHY", FiniteStrainBehaviourOptions::CAUCHY);
  m.set_const("PK1", FiniteStrainBehaviourOptions::PK1);
  m.set_const("FirstPiolaKirchhoffStress", FiniteStrainBehaviourOptions::PK1);
  m.set_const("PK2", FiniteStrainBehaviourOptions::PK2);
  m.set_const("SecondPiolaKirchhoffStress", FiniteStrainBehaviourOptions::PK2);
  //!
  m.add_bits<FiniteStrainBehaviourOptions::TangentOperator>(
      "FiniteStrainBehaviourOptionsTangentOperator");
  m.set_const("DSIG_DF", FiniteStrainBehaviourOptions::DSIG_DF);
  m.set_const("DCAUCHY_DF", FiniteStrainBehaviourOptions::DSIG_DF);
  m.set_const("DPK1_DF", FiniteStrainBehaviourOptions::DPK1_DF);
  m.set_const("DS_DEGL", FiniteStrainBehaviourOptions::DS_DEGL);
  m.set_const("DPK2_DEGL", FiniteStrainBehaviourOptions::DS_DEGL);
  //!
  m.add_bits<Behaviour::Symmetry>("BehaviourSymmetry");
  m.set_const("ISOTROPIC", Behaviour::ISOTROPIC);
  m.set_const("Isotropic", Behaviour::ISOTROPIC);
  m.set_const("ORTHOTROPIC", Behaviour::ORTHOTROPIC);
  m.set_const("Orthotropic", Behaviour::ORTHOTROPIC);
  //!
  m.add_type<FiniteStrainBehaviourOptions>("FiniteStrainBehaviourOptions")
      .method("get_stress_measure",
              [](const FiniteStrainBehaviourOptions &o) {
                return o.stress_measure;
              })
      .method("set_stress_measure!",
              [](FiniteStrainBehaviourOptions &o,
                 const FiniteStrainBehaviourOptions::StressMeasure &s) {
                o.stress_measure = s;
              })
      .method("get_tangent_operator",
              [](const FiniteStrainBehaviourOptions &o) {
                return o.tangent_operator;
              })
      .method("set_tangent_operator!",
              [](FiniteStrainBehaviourOptions &o,
                 const FiniteStrainBehaviourOptions::TangentOperator &to) {
                o.tangent_operator = to;
              });
  //
  m.add_bits<Behaviour::BehaviourType>("BehaviourType");
  m.set_const("GENERALBEHAVIOUR", Behaviour::GENERALBEHAVIOUR);
  m.set_const("GeneralBehaviour", Behaviour::GENERALBEHAVIOUR);
  m.set_const("STANDARDSTRAINBASEDBEHAVIOUR",
              Behaviour::STANDARDSTRAINBASEDBEHAVIOUR);
  m.set_const("StandardStrainBasedBehaviour",
              Behaviour::STANDARDSTRAINBASEDBEHAVIOUR);
  m.set_const("STANDARDFINITESTRAINBEHAVIOUR",
              Behaviour::STANDARDFINITESTRAINBEHAVIOUR);
  m.set_const("StandardFiniteStrainBehaviour",
              Behaviour::STANDARDFINITESTRAINBEHAVIOUR);
  m.set_const("COHESIVEZONEMODEL", Behaviour::COHESIVEZONEMODEL);
  m.set_const("CohesiveZoneModel", Behaviour::COHESIVEZONEMODEL);
  //! kinematic of the behaviour treated
  m.add_bits<Behaviour::Kinematic>("BehaviourKinematic");
  m.set_const("UNDEFINEDKINEMATIC", Behaviour::UNDEFINEDKINEMATIC);
  m.set_const("UndefinedKinematic", Behaviour::UNDEFINEDKINEMATIC);
  m.set_const("SMALLSTRAINKINEMATIC", Behaviour::SMALLSTRAINKINEMATIC);
  m.set_const("SmallStrainKinematic", Behaviour::SMALLSTRAINKINEMATIC);
  m.set_const("COHESIVEZONEKINEMATIC", Behaviour::COHESIVEZONEKINEMATIC);
  m.set_const("CohesiveZoneKinematic", Behaviour::COHESIVEZONEKINEMATIC);
  m.set_const("FINITESTRAINKINEMATIC_ETO_PK1",
              Behaviour::FINITESTRAINKINEMATIC_ETO_PK1);
  m.set_const("FiniteStrainKinematic_ETO_PK1",
              Behaviour::FINITESTRAINKINEMATIC_ETO_PK1);
  m.set_const("FINITESTRAINKINEMATIC_F_CAUCHY",
              Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY);
  m.set_const("FiniteStrainKinematic_F_CAUCHY",
              Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY);
  //
  m.add_type<Behaviour>("Behaviour")
      .method("get_library", [](const Behaviour &b) { return b.library; })
      .method("get_behaviour", [](const Behaviour &b) { return b.behaviour; })
      .method("get_hypothesis", [](const Behaviour &b) { return b.hypothesis; })
      .method("get_function", [](const Behaviour &b) { return b.function; })
      .method("get_source", [](const Behaviour &b) { return b.source; })
      .method("get_tfel_version",
              [](const Behaviour &b) { return b.tfel_version; })
      .method("get_behaviour_type", [](const Behaviour &b) { return b.btype; })
      .method("get_symmetry", [](const Behaviour &b) { return b.symmetry; })
      .method("get_kinematic", [](const Behaviour &b) { return b.kinematic; })
      .method("get_material_properties",
              [](const Behaviour &b) { return b.mps; })
      .method("get_internal_state_variables",
              [](const Behaviour &b) { return b.isvs; })
      .method("get_external_state_variables",
              [](const Behaviour &b) { return b.esvs; })
      .method("get_parameters", [](const Behaviour &b) { return b.params; });
  //
  m.method("load", load);
  m.method("load", load2);
  m.method("set_parameter!", setParameter1);
  m.method("set_integer_parameter!", setParameter2);
  m.method("set_unsigned_short_parameter!",
           [](const Behaviour &b, const std::string &n, const unsigned int v) {
             mgis::behaviour::setParameter(b, n,
                                           static_cast<unsigned short>(v));
           });
  m.method("get_parameter_default_value", getParameterDefaultValue1);
  m.method("get_integer_parameter_default_value", getParameterDefaultValue2);
  m.method("get_unsigned_short_parameter_default_value",
           getParameterDefaultValue3);
  m.method("has_bounds", hasBounds);
  m.method("has_lower_bound", hasLowerBound);
  m.method("has_upper_bound", hasUpperBound);
  m.method("get_lower_bound", getLowerBound);
  m.method("get_upper_bound", getUpperBound);
  m.method("has_physical_bounds", hasPhysicalBounds);
  m.method("has_lower_physical_bound", hasLowerPhysicalBound);
  m.method("has_upper_physical_bound", hasUpperPhysicalBound);
  m.method("get_lower_physical_bound", getLowerPhysicalBound);
  m.method("get_upper_physical_bound", getUpperPhysicalBound);

}  // end of declareBehaviour
