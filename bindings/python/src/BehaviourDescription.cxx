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

#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/def.hpp>
#include "MGIS/Python/VectorConverter.hxx"
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"

// forward declaration
void declareBehaviourDescription();

static const char *BehaviourDescription_getType(
    const mgis::behaviour::BehaviourDescription &b) {
  using mgis::behaviour::BehaviourDescription;
  switch (b.btype) {
    case BehaviourDescription::GENERALBEHAVIOUR:
      return "GeneralBehaviour";
      break;
    case BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR:
      return "StandardStrainBasedBehaviour";
      break;
    case BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR:
      return "StandardFiniteStrainBehaviour";
      break;
    case BehaviourDescription::COHESIVEZONEMODEL:
      return "CohesiveZoneModel";
      break;
    default:
      mgis::raise("BehaviourDescription_getType: unsupported behaviour type");
  }
  return "";
}  // end of BehaviourDescription_getType

static const char *BehaviourDescription_getKinematic(
    const mgis::behaviour::BehaviourDescription &b) {
  using mgis::behaviour::BehaviourDescription;
  switch (b.kinematic) {
    case BehaviourDescription::SMALLSTRAINKINEMATIC:
      return "SmallStrainKinematic";
      break;
    case BehaviourDescription::COHESIVEZONEKINEMATIC:
      return "CohesiveZoneKinematic";
      break;
    case BehaviourDescription::FINITESTRAINKINEMATIC_F_CAUCHY:
      return "F_CAUCHY";
      break;
    case BehaviourDescription::FINITESTRAINKINEMATIC_ETO_PK1:
      return "ETO_PK1";
      break;
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
      mgis::raise("unsupported symmetry type");
      break;
  }
  return "UndefinedSymemtry";
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

static boost::python::list BehaviourDescription_getParameters(
    const mgis::behaviour::BehaviourDescription &b) {
  return mgis::python::convert_vector_to_list(b.params);
}  // end of BehaviourDescription_getParameters

static boost::python::list BehaviourDescription_getIntegerParameters(
    const mgis::behaviour::BehaviourDescription &b) {
  return mgis::python::convert_vector_to_list(b.iparams);
}  // end of BehaviourDescription_getIntegerParameters

static boost::python::list BehaviourDescription_getUnsignedShortParameters(
    const mgis::behaviour::BehaviourDescription &b) {
  return mgis::python::convert_vector_to_list(b.usparams);
}  // end of BehaviourDescription_getUnsignedShortParameters

static std::vector<
    std::pair<mgis::behaviour::Variable, mgis::behaviour::Variable>>
BehaviourDescription_getTangentOperatorBlocks(
    const mgis::behaviour::BehaviourDescription &b) {
  return b.to_blocks;
}  // end of BehaviourDescription_getTangentOperatorBlocks

void declareBehaviourDescription() {
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
  boost::python::enum_<BehaviourDescription::Symmetry>("BehaviourSymmetry")
      .value("ISOTROPIC", BehaviourDescription::Symmetry::ISOTROPIC)
      .value("Isotropic", BehaviourDescription::Symmetry::ISOTROPIC)
      .value("ORTHOTROPIC", BehaviourDescription::Symmetry::ORTHOTROPIC)
      .value("Orthotropic", BehaviourDescription::Symmetry::ORTHOTROPIC);
  // wrapping the BehaviourDescription::BehaviourType enum
  boost::python::enum_<BehaviourDescription::BehaviourType>("BehaviourType")
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
  boost::python::enum_<BehaviourDescription::Kinematic>("BehaviourKinematic")
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
  boost::python::class_<BehaviourDescription, boost::noncopyable>(
      "BehaviourDescription", boost::python::no_init)
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
      .def_readonly("tfel_version", &BehaviourDescription::tfel_version,
                    "version of TFEL used to generate the behaviour")
      .add_property("btype", &BehaviourDescription::btype,
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
      .add_property("kinematic", &BehaviourDescription::kinematic,
                    "return the behaviour kinematic")
      .def("getKinematic", &BehaviourDescription_getKinematic,
           "return the behaviour kinematic")
      .def("getSymmetry", &BehaviourDescription_getSymmetry,
           "return the behaviour symmetry")
      .add_property("symmetry", &BehaviourDescription::symmetry,
                    "return the behaviour symmetry")
      .add_property("gradients", &BehaviourDescription_getGradients,
                    "list of gradients")
      .add_property("thermodynamic_forces",
                    &BehaviourDescription_getThermodynamicForces,
                    "list of thermodynamic forces")
      .add_property("mps", &BehaviourDescription_getMaterialProperties,
                    "list of material properties")
      .add_property("material_properties",
                    &BehaviourDescription_getMaterialProperties,
                    "list of material properties (same as the `mps` property)")
      .add_property("isvs", &BehaviourDescription_getInternalStateVariables,
                    "list of internal state variables")
      .add_property(
          "internal_state_variables",
          &BehaviourDescription_getInternalStateVariables,
          "list of internal state variables (same as the `isvs` property)")
      .add_property("esvs", &BehaviourDescription_getExternalStateVariables,
                    "list of external state variables")
      .add_property(
          "external_state_variables",
          &BehaviourDescription_getExternalStateVariables,
          "list of external state variables (same as the `esvs` property)")
      .add_property("params", &BehaviourDescription_getParameters,
                    "list of parameters")
      .add_property("parameters", &BehaviourDescription_getParameters,
                    "list of parameters (same as the `params` property")
      .add_property("iparams", &BehaviourDescription_getIntegerParameters,
                    "list of integer parameters")
      .add_property(
          "integer_parameters", &BehaviourDescription_getIntegerParameters,
          "list of integer parameters (same as the `integer_params` property")
      .add_property("usparams",
                    &BehaviourDescription_getUnsignedShortParameters,
                    "list of unsigned short parameters")
      .add_property(
          "unsigned_short_parameters",
          &BehaviourDescription_getUnsignedShortParameters,
          "list of unsigned short parameters (same as the `usparams` property)")
      .add_property("tangent_operator_blocks",
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

  boost::python::def(
      "isStandardFiniteStrainBehaviour",
      mgis::behaviour::isStandardFiniteStrainBehaviour,
      "return if the given behaviour is a standard finite strain behaviour, "
      "i.e. is a finite strain behaviour using the standard finite strain "
      "kinematic (called F-Cauchy although the stress measure can be chosen "
      "when loading the behaviour)");
  boost::python::def("getParameterDefaultValue", getParameterDefaultValue1);
  boost::python::def("getIntegerParameterDefaultValue",
                     getParameterDefaultValue2);
  boost::python::def("getUnsignedShortParameterDefaultValue",
                     getParameterDefaultValue3);
  boost::python::def("hasBounds", hasBounds);
  boost::python::def("hasLowerBound", hasLowerBound);
  boost::python::def("hasUpperBound", hasUpperBound);
  boost::python::def("getLowerBound", getLowerBound);
  boost::python::def("getUpperBound", getUpperBound);
  boost::python::def("hasPhysicalBounds", hasPhysicalBounds);
  boost::python::def("hasLowerPhysicalBound", hasLowerPhysicalBound);
  boost::python::def("hasUpperPhysicalBound", hasUpperPhysicalBound);
  boost::python::def("getLowerPhysicalBound", getLowerPhysicalBound);
  boost::python::def("getUpperPhysicalBound", getUpperPhysicalBound);

}  // end of declareBehaviour
