/*!
 * \file   bindings/python/src/Behaviour.cxx
 * \brief
 * \author Thomas Helfer
 * \date   06/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/def.hpp>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

// forward declaration
void declareBehaviour();

static const char* Behaviour_getType(const mgis::behaviour::Behaviour& b) {
  using mgis::behaviour::Behaviour;
  switch (b.btype) {
    case Behaviour::GENERALBEHAVIOUR:
      return "GeneralBehaviour";
      break;
    case Behaviour::STANDARDSTRAINBASEDBEHAVIOUR:
      return "StandardStrainBasedBehaviour";
      break;
    case Behaviour::STANDARDFINITESTRAINBEHAVIOUR:
      return "StandardFiniteStrainBehaviour";
      break;
    case Behaviour::COHESIVEZONEMODEL:
      return "CohesiveZoneModel";
      break;
    default:
      mgis::raise("Behaviour_getType: unsupported behaviour type");
  }
  return "";
}  // end of Behaviour_getType

static const char* Behaviour_getKinematic(const mgis::behaviour::Behaviour& b) {
  using mgis::behaviour::Behaviour;
  switch (b.kinematic) {
    case Behaviour::SMALLSTRAINKINEMATIC:
      return "SmallStrainKinematic";
      break;
    case Behaviour::COHESIVEZONEKINEMATIC:
      return "CohesiveZoneKinematic";
      break;
    case Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY:
      return "F_CAUCHY";
      break;
    case Behaviour::FINITESTRAINKINEMATIC_ETO_PK1:
      return "ETO_PK1";
      break;
    case Behaviour::UNDEFINEDKINEMATIC:
    default:
      break;
  }
  return "UndefinedKinematic";
}  // end of Behaviour_getKinematic

static const char* Behaviour_getSymmetry(const mgis::behaviour::Behaviour& b) {
  using mgis::behaviour::Behaviour;
  switch (b.kinematic) {
    case Behaviour::ISOTROPIC:
      return "Isotropic";
      break;
    case Behaviour::ORTHOTROPIC:
      return "Orthotropic";
      break;
    default:
      mgis::raise("unsupported symmetry type");
      break;
  }
  return "UndefinedSymemtry";
}  // end of Behaviour_getSymmetry

static std::vector<mgis::behaviour::Variable> Behaviour_getGradients(
    const mgis::behaviour::Behaviour& b) {
  return b.gradients;
}  // end of Behaviour_getGradients

static std::vector<mgis::behaviour::Variable> Behaviour_getThermodynamicForces(
    const mgis::behaviour::Behaviour& b) {
  return b.thermodynamic_forces;
}  // end of Behaviour_getThermodynamicForces

static std::vector<mgis::behaviour::Variable> Behaviour_getMaterialProperties(
    const mgis::behaviour::Behaviour& b) {
  return b.mps;
}  // end of Behaviour_getMaterialProperties

static std::vector<mgis::behaviour::Variable>
Behaviour_getInternalStateVariables(const mgis::behaviour::Behaviour& b) {
  return b.isvs;
}  // end of Behaviour_getInternalStateVariables

static std::vector<mgis::behaviour::Variable>
Behaviour_getExternalStateVariables(const mgis::behaviour::Behaviour& b) {
  return b.esvs;
}  // end of Behaviour_getExternalStateVariables

static std::vector<mgis::behaviour::Variable> Behaviour_getParameters(
    const mgis::behaviour::Behaviour& b) {
  return b.parameters;
}  // end of Behaviour_getParameters

void declareBehaviour() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::Hypothesis;
  Behaviour (*load_ptr)(const std::string&, const std::string&,
                        const Hypothesis) = &mgis::behaviour::load;
  // wrapping the Behaviour::Symmetry enum
  boost::python::enum_<Behaviour::Symmetry>("BehaviourSymmetry")
      .value("ISOTROPIC", Behaviour::Symmetry::ISOTROPIC)
      .value("Isotropic", Behaviour::Symmetry::ISOTROPIC)
      .value("ORTHOTROPIC", Behaviour::Symmetry::ORTHOTROPIC)
      .value("Orthotropic", Behaviour::Symmetry::ORTHOTROPIC);
  // wrapping the Behaviour::BehaviourType enum
  boost::python::enum_<Behaviour::BehaviourType>("BehaviourType")
      .value("GENERALBEHAVIOUR", Behaviour::BehaviourType::GENERALBEHAVIOUR)
      .value("GeneralBehaviour", Behaviour::BehaviourType::GENERALBEHAVIOUR)
      .value("STANDARDSTRAINBASEDBEHAVIOUR",
             Behaviour::BehaviourType::STANDARDSTRAINBASEDBEHAVIOUR)
      .value("StandardStrainBasedBehaviour",
             Behaviour::BehaviourType::STANDARDSTRAINBASEDBEHAVIOUR)
      .value("STANDARDFINITESTRAINBEHAVIOUR",
             Behaviour::BehaviourType::STANDARDFINITESTRAINBEHAVIOUR)
      .value("StandardFiniteStrainBehaviour",
             Behaviour::BehaviourType::STANDARDFINITESTRAINBEHAVIOUR)
      .value("COHESIVEZONEMODEL", Behaviour::BehaviourType::COHESIVEZONEMODEL)
      .value("CohesiveZoneModel", Behaviour::BehaviourType::COHESIVEZONEMODEL);
  // wrapping the Behaviour::BehaviourType enum
  boost::python::enum_<Behaviour::Kinematic>("BehaviourKinematic")
      .value("UNDEFINEDKINEMATIC", Behaviour::Kinematic::UNDEFINEDKINEMATIC)
      .value("UndefinedKinematic", Behaviour::Kinematic::UNDEFINEDKINEMATIC)
      .value("SMALLSTRAINKINEMATIC", Behaviour::Kinematic::SMALLSTRAINKINEMATIC)
      .value("SmallStrainKinematic", Behaviour::Kinematic::SMALLSTRAINKINEMATIC)
      .value("COHESIVEZONEKINEMATIC",
             Behaviour::Kinematic::COHESIVEZONEKINEMATIC)
      .value("CohesiveZoneKinematic",
             Behaviour::Kinematic::COHESIVEZONEKINEMATIC)
      .value("FINITESTRAINKINEMATIC_F_CAUCHY",
             Behaviour::Kinematic::FINITESTRAINKINEMATIC_F_CAUCHY)
      .value("FiniteStrainKinematic_F_Cauchy",
             Behaviour::Kinematic::FINITESTRAINKINEMATIC_F_CAUCHY)
      .value("FINITESTRAINKINEMATIC_Eto_PK1",
             Behaviour::Kinematic::FINITESTRAINKINEMATIC_ETO_PK1)
      .value("FiniteStrainKinematic_Eto_PK1",
             Behaviour::Kinematic::FINITESTRAINKINEMATIC_ETO_PK1);
  // wrapping the Behaviour class
  boost::python::class_<Behaviour>("Behaviour")
      .def_readonly("library", &Behaviour::library,
                    "name of the library in which the behaviour is implemented")
      .def_readonly("behaviour", &Behaviour::behaviour, "name of the behaviour")
      .def_readonly("hypothesis", &Behaviour::hypothesis,
                    "modelling hypothesis")
      .def_readonly("function", &Behaviour::function,
                    "function implementing the behaviour")
      .def_readonly("source", &Behaviour::source,
                    "name of the `MFront` source file")
      .def_readonly("tfel_version", &Behaviour::tfel_version,
                    "version of TFEL used to generate the behaviour")
      .add_property("btype", &Behaviour::btype,
                    "return the type of the behaviour")
      .def("getBehaviourType", &Behaviour_getType,
           "return the type of the behaviour")
      .add_property("kinematic", &Behaviour::kinematic,
                    "return the behaviour kinematic")
      .def("getKinematic", &Behaviour_getKinematic,
           "return the behaviour kinematic")
      .def("getSymmetry", &Behaviour_getSymmetry, "return the behaviour symmetry")
      .add_property("symmetry", &Behaviour::symmetry,
                    "return the behaviour symmetry")
      .add_property("gradients", &Behaviour_getGradients, "list of gradients")
      .add_property("thermodynamic_forces", &Behaviour_getThermodynamicForces,
                    "list of thermodynamic forces")
      .add_property("mps", &Behaviour_getMaterialProperties,
                    "list of material properties")
      .add_property("material_properties", &Behaviour_getMaterialProperties,
                    "list of material properties (same as the `mps` property)")
      .add_property("isvs", &Behaviour_getInternalStateVariables,
                    "list of internal state variables")
      .add_property(
          "internal_state_variables", &Behaviour_getInternalStateVariables,
          "list of internal state variables (same as the `isvs` property)")
      .add_property("esvs", &Behaviour_getExternalStateVariables,
                    "list of external state variables")
      .add_property(
          "external_state_variables", &Behaviour_getExternalStateVariables,
          "list of external state variables (same as the `esvs` property)")
      .add_property("parameters", &Behaviour_getParameters,
                    "list of parameters");

  // wrapping the load function
  boost::python::def("load", load_ptr);
}  // end of declareBehaviour
