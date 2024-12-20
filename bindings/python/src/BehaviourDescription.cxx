/*!
 * \file   bindings/python/src/BehaviourDescription.cxx
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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"

// forward declaration
void declareBehaviourDescription(pybind11::module_&);

static const char *BehaviourDescription_getType(
    const mgis::behaviour::BehaviourDescription &b) {
  using mgis::behaviour::BehaviourDescription;
  switch (b.btype) {
    case BehaviourDescription::GENERALBEHAVIOUR:
      return "GeneralBehaviour";
    case BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR:
      return "StandardStrainBasedBehaviour";
    case BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR:
      return "StandardFiniteStrainBehaviour";
    case BehaviourDescription::COHESIVEZONEMODEL:
      return "CohesiveZoneModel";
    default:
      break;
  }
  mgis::raise("BehaviourDescription_getType: unsupported behaviour type");
}  // end of BehaviourDescription_getType

static const char *BehaviourDescription_getKinematic(
    const mgis::behaviour::BehaviourDescription &b) {
  using mgis::behaviour::BehaviourDescription;
  switch (b.kinematic) {
    case BehaviourDescription::SMALLSTRAINKINEMATIC:
      return "SmallStrainKinematic";
    case BehaviourDescription::COHESIVEZONEKINEMATIC:
      return "CohesiveZoneKinematic";
    case BehaviourDescription::FINITESTRAINKINEMATIC_F_CAUCHY:
      return "F_CAUCHY";
    case BehaviourDescription::FINITESTRAINKINEMATIC_ETO_PK1:
      return "ETO_PK1";
    case BehaviourDescription::UNDEFINEDKINEMATIC:
    default:
      break;
  }
  return "UndefinedKinematic";
}  // end of BehaviourDescription_getKinematic

static const char *BehaviourDescription_getSymmetry(
    const mgis::behaviour::BehaviourDescription &b) {
  using mgis::behaviour::BehaviourDescription;
  switch (b.symmetry) {
    case BehaviourDescription::ISOTROPIC:
      return "Isotropic";
    case BehaviourDescription::ORTHOTROPIC:
      return "Orthotropic";
    default:
      break;
  }
  mgis::raise("unsupported symmetry type");
}  // end of BehaviourDescription_getSymmetry

static std::vector<mgis::behaviour::Variable> BehaviourDescription_getGradients(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.gradients;
}  // end of BehaviourDescription_getGradients

static std::vector<mgis::behaviour::Variable>
BehaviourDescription_getThermodynamicForces(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.thermodynamic_forces;
}  // end of BehaviourDescription_getThermodynamicForces

static std::vector<mgis::behaviour::Variable>
BehaviourDescription_getMaterialProperties(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.mps;
}  // end of BehaviourDescription_getMaterialProperties

static std::vector<mgis::behaviour::Variable>
BehaviourDescription_getInternalStateVariables(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.isvs;
}  // end of BehaviourDescription_getInternalStateVariables

static std::vector<mgis::behaviour::Variable>
BehaviourDescription_getExternalStateVariables(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.esvs;
}  // end of BehaviourDescription_getExternalStateVariables

static std::vector<
    std::pair<mgis::behaviour::Variable, mgis::behaviour::Variable>>
BehaviourDescription_getTangentOperatorBlocks(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.to_blocks;
}  // end of BehaviourDescription_getTangentOperatorBlocks

void declareBehaviourDescription(pybind11::module_ &m) {
  using mgis::behaviour::BehaviourDescription;
  //
  double (*getParameterDefaultValue1)(const BehaviourDescription &,
                                      const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<double>;
  int (*getParameterDefaultValue2)(const BehaviourDescription &,
                                   const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<int>;
  unsigned short (*getParameterDefaultValue3)(const BehaviourDescription &,
                                              const std::string &) =
      &mgis::behaviour::getParameterDefaultValue<unsigned short>;
  bool (*hasBounds)(const BehaviourDescription &, const std::string &) =
      mgis::behaviour::hasBounds;
  bool (*hasLowerBound)(const BehaviourDescription &, const std::string &) =
      mgis::behaviour::hasLowerBound;
  bool (*hasUpperBound)(const BehaviourDescription &, const std::string &) =
      mgis::behaviour::hasUpperBound;
  long double (*getLowerBound)(const BehaviourDescription &,
                               const std::string &) =
      mgis::behaviour::getLowerBound;
  long double (*getUpperBound)(const BehaviourDescription &,
                               const std::string &) =
      mgis::behaviour::getUpperBound;
  bool (*hasPhysicalBounds)(const BehaviourDescription &, const std::string &) =
      mgis::behaviour::hasPhysicalBounds;
  bool (*hasLowerPhysicalBound)(const BehaviourDescription &,
                                const std::string &) =
      mgis::behaviour::hasLowerPhysicalBound;
  bool (*hasUpperPhysicalBound)(const BehaviourDescription &,
                                const std::string &) =
      mgis::behaviour::hasUpperPhysicalBound;
  long double (*getLowerPhysicalBound)(const BehaviourDescription &,
                                       const std::string &) =
      mgis::behaviour::getLowerPhysicalBound;
  long double (*getUpperPhysicalBound)(const BehaviourDescription &,
                                       const std::string &) =
      mgis::behaviour::getUpperPhysicalBound;

  // wrapping the BehaviourDescription::Symmetry enum
  pybind11::enum_<BehaviourDescription::Symmetry>(m, "BehaviourSymmetry")
      .value("ISOTROPIC", BehaviourDescription::Symmetry::ISOTROPIC)
      .value("Isotropic", BehaviourDescription::Symmetry::ISOTROPIC)
      .value("ORTHOTROPIC", BehaviourDescription::Symmetry::ORTHOTROPIC)
      .value("Orthotropic", BehaviourDescription::Symmetry::ORTHOTROPIC);
  // wrapping the BehaviourDescription::BehaviourType enum
  pybind11::enum_<BehaviourDescription::BehaviourType>(m, "BehaviourType")
      .value("GENERALBEHAVIOUR",
             BehaviourDescription::BehaviourType::GENERALBEHAVIOUR)
      .value("GeneralBehaviour",
             BehaviourDescription::BehaviourType::GENERALBEHAVIOUR)
      .value("STANDARDSTRAINBASEDBEHAVIOUR",
             BehaviourDescription::BehaviourType::STANDARDSTRAINBASEDBEHAVIOUR)
      .value("StandardStrainBasedBehaviour",
             BehaviourDescription::BehaviourType::STANDARDSTRAINBASEDBEHAVIOUR)
      .value("STANDARDFINITESTRAINBEHAVIOUR",
             BehaviourDescription::BehaviourType::STANDARDFINITESTRAINBEHAVIOUR)
      .value("StandardFiniteStrainBehaviour",
             BehaviourDescription::BehaviourType::STANDARDFINITESTRAINBEHAVIOUR)
      .value("COHESIVEZONEMODEL",
             BehaviourDescription::BehaviourType::COHESIVEZONEMODEL)
      .value("CohesiveZoneModel",
             BehaviourDescription::BehaviourType::COHESIVEZONEMODEL);
  // wrapping the BehaviourDescription::BehaviourType enum
  pybind11::enum_<BehaviourDescription::Kinematic>(m, "BehaviourKinematic")
      .value("UNDEFINEDKINEMATIC",
             BehaviourDescription::Kinematic::UNDEFINEDKINEMATIC)
      .value("UndefinedKinematic",
             BehaviourDescription::Kinematic::UNDEFINEDKINEMATIC)
      .value("SMALLSTRAINKINEMATIC",
             BehaviourDescription::Kinematic::SMALLSTRAINKINEMATIC)
      .value("SmallStrainKinematic",
             BehaviourDescription::Kinematic::SMALLSTRAINKINEMATIC)
      .value("COHESIVEZONEKINEMATIC",
             BehaviourDescription::Kinematic::COHESIVEZONEKINEMATIC)
      .value("CohesiveZoneKinematic",
             BehaviourDescription::Kinematic::COHESIVEZONEKINEMATIC)
      .value("FINITESTRAINKINEMATIC_F_CAUCHY",
             BehaviourDescription::Kinematic::FINITESTRAINKINEMATIC_F_CAUCHY)
      .value("FiniteStrainKinematic_F_Cauchy",
             BehaviourDescription::Kinematic::FINITESTRAINKINEMATIC_F_CAUCHY)
      .value("FINITESTRAINKINEMATIC_Eto_PK1",
             BehaviourDescription::Kinematic::FINITESTRAINKINEMATIC_ETO_PK1)
      .value("FiniteStrainKinematic_Eto_PK1",
             BehaviourDescription::Kinematic::FINITESTRAINKINEMATIC_ETO_PK1);
  // wrapping the Behaviour class
  pybind11::class_<BehaviourDescription>(m, "BehaviourDescription")
      .def_readonly("library", &BehaviourDescription::library,
                    "name of the library in which the behaviour is implemented")
      .def_readonly("behaviour", &BehaviourDescription::behaviour,
                    "name of the behaviour")
      .def_readonly("hypothesis", &BehaviourDescription::hypothesis,
                    "modelling hypothesis")
      .def_readonly("function", &BehaviourDescription::function,
                    "function implementing the behaviour")
      .def_readonly("source", &BehaviourDescription::source,
                    "name of the `MFront` source file")
      .def_readonly("author", &BehaviourDescription::author,
                    "author of the `MFront`'s file")
      .def_readonly("date", &BehaviourDescription::date, "date")
      .def_readonly("validator", &BehaviourDescription::validator,
                    "validator of the `MFront`'s file")
      .def_readonly("build_id", &BehaviourDescription::build_id,
                    "build identifier of the `MFront`'s file")
      .def_readonly("build_identifier", &BehaviourDescription::build_id,
                    "build identifier of the `MFront`'s file")
      .def_readonly("tfel_version", &BehaviourDescription::tfel_version,
                    "version of TFEL used to generate the behaviour")
      .def_readonly("btype", &BehaviourDescription::btype,
                    "return the type of the behaviour")
      .def_readonly(
          "computesStoredEnergy", &BehaviourDescription::computesStoredEnergy,
          "a boolean stating if the behaviour computes the stored energy")
      .def_readonly(
          "computesDissipatedEnergy",
          &BehaviourDescription::computesDissipatedEnergy,
          "a boolean stating if the behaviour computes the dissipated energy")
      .def("getBehaviourType", &BehaviourDescription_getType,
           "return the type of the behaviour")
      .def_readonly("kinematic", &BehaviourDescription::kinematic,
                    "return the behaviour kinematic")
      .def("getKinematic", &BehaviourDescription_getKinematic,
           "return the behaviour kinematic")
      .def("getSymmetry", &BehaviourDescription_getSymmetry,
           "return the behaviour symmetry")
      .def_readonly("symmetry", &BehaviourDescription::symmetry,
                    "return the behaviour symmetry")
      .def_property_readonly("gradients", &BehaviourDescription_getGradients,
                             "list of gradients")
      .def_property_readonly("thermodynamic_forces",
                             &BehaviourDescription_getThermodynamicForces,
                             "list of thermodynamic forces")
      .def_property_readonly("mps", &BehaviourDescription_getMaterialProperties,
                             "list of material properties")
      .def_property_readonly(
          "material_properties", &BehaviourDescription_getMaterialProperties,
          "list of material properties (same as the `mps` property)")
      .def_property_readonly("isvs",
                             &BehaviourDescription_getInternalStateVariables,
                             "list of internal state variables")
      .def_property_readonly(
          "internal_state_variables",
          &BehaviourDescription_getInternalStateVariables,
          "list of internal state variables (same as the `isvs` property)")
      .def_property_readonly("esvs",
                             &BehaviourDescription_getExternalStateVariables,
                             "list of external state variables")
      .def_property_readonly(
          "external_state_variables",
          &BehaviourDescription_getExternalStateVariables,
          "list of external state variables (same as the `esvs` property)")
      .def_readonly("params", &BehaviourDescription::params,
                    "list of parameters")
      .def_readonly("parameters", &BehaviourDescription::params,
                    "list of parameters (same as the `params` property")
      .def_readonly("iparams", &BehaviourDescription::iparams,
                    "list of integer parameters")
      .def_readonly(
          "integer_parameters", &BehaviourDescription::iparams,
          "list of integer parameters (same as the `integer_params` property")
      .def_readonly("usparams", &BehaviourDescription::usparams,
                    "list of unsigned short parameters")
      .def_readonly(
          "unsigned_short_parameters", &BehaviourDescription::usparams,
          "list of unsigned short parameters (same as the `usparams` property)")
      .def_property_readonly("tangent_operator_blocks",
                             BehaviourDescription_getTangentOperatorBlocks)
      .def("getParameterDefaultValue", getParameterDefaultValue1)
      .def("getIntegerParameterDefaultValue", getParameterDefaultValue2)
      .def("getUnsignedShortParameterDefaultValue", getParameterDefaultValue3)
      .def("hasBounds", hasBounds)
      .def("hasLowerBound", hasLowerBound)
      .def("hasUpperBound", hasUpperBound)
      .def("getLowerBound ", getLowerBound)
      .def("getUpperBound", getUpperBound)
      .def("hasPhysicalBounds", hasPhysicalBounds)
      .def("hasLowerPhysicalBound", hasLowerPhysicalBound)
      .def("hasUpperPhysicalBound", hasUpperPhysicalBound)
      .def("getLowerPhysicalBound ", getLowerPhysicalBound)
      .def("getUpperPhysicalBound ", getUpperPhysicalBound);

  m.def(
      "isStandardFiniteStrainBehaviour",
      mgis::behaviour::isStandardFiniteStrainBehaviour,
      "return if the given behaviour is a standard finite strain behaviour, "
      "i.e. is a finite strain behaviour using the standard finite strain "
      "kinematic (called F-Cauchy although the stress measure can be chosen "
      "when loading the behaviour)");
  m.def("getParameterDefaultValue", getParameterDefaultValue1);
  m.def("getIntegerParameterDefaultValue",
                     getParameterDefaultValue2);
  m.def("getUnsignedShortParameterDefaultValue",
                     getParameterDefaultValue3);
  m.def("hasBounds", hasBounds);
  m.def("hasLowerBound", hasLowerBound);
  m.def("hasUpperBound", hasUpperBound);
  m.def("getLowerBound", getLowerBound);
  m.def("getUpperBound", getUpperBound);
  m.def("hasPhysicalBounds", hasPhysicalBounds);
  m.def("hasLowerPhysicalBound", hasLowerPhysicalBound);
  m.def("hasUpperPhysicalBound", hasUpperPhysicalBound);
  m.def("getLowerPhysicalBound", getLowerPhysicalBound);
  m.def("getUpperPhysicalBound", getUpperPhysicalBound);

}  // end of declareBehaviour
