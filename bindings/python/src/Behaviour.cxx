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
#include "MGIS/Python/VectorConverter.hxx"
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
  switch (b.symmetry) {
    case Behaviour::ISOTROPIC:
      return "Isotropic";
    case Behaviour::ORTHOTROPIC:
      return "Orthotropic";
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

static boost::python::list Behaviour_getParameters(
    const mgis::behaviour::Behaviour& b) {
  return mgis::python::convert_vector_to_list(b.params);
}  // end of Behaviour_getParameters

static boost::python::list Behaviour_getIntegerParameters(
    const mgis::behaviour::Behaviour& b) {
  return mgis::python::convert_vector_to_list(b.iparams);
}  // end of Behaviour_getIntegerParameters

static boost::python::list Behaviour_getUnsignedShortParameters(
    const mgis::behaviour::Behaviour& b) {
  return mgis::python::convert_vector_to_list(b.usparams);
}  // end of Behaviour_getUnsignedShortParameters

static std::vector<
    std::pair<mgis::behaviour::Variable, mgis::behaviour::Variable>>
Behaviour_getTangentOperatorBlocks(const mgis::behaviour::Behaviour &b) {
  return b.to_blocks;
}  // end of Behaviour_getUnsignedShortParameters


void declareBehaviour() {
  using mgis::behaviour::FiniteStrainBehaviourOptions;
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::Hypothesis;
  Behaviour (*load_ptr)(const std::string&, const std::string&,
                        const Hypothesis) = &mgis::behaviour::load;
  Behaviour (*load_ptr2)(const FiniteStrainBehaviourOptions&,
			 const std::string&, const std::string&,
			 const Hypothesis) = &mgis::behaviour::load;
  void (*setParameter1)(const Behaviour&, const std::string&, const double) =
      &mgis::behaviour::setParameter;
  void (*setParameter2)(const Behaviour&, const std::string&, const int) =
      &mgis::behaviour::setParameter;
  void (*setParameter3)(const Behaviour&, const std::string&, const unsigned short) =
      &mgis::behaviour::setParameter;
  double (*getParameterDefaultValue1)(const Behaviour&, const std::string&) =
      &mgis::behaviour::getParameterDefaultValue<double>;
  int (*getParameterDefaultValue2)(const Behaviour&, const std::string&) =
      &mgis::behaviour::getParameterDefaultValue<int>;
  unsigned short (*getParameterDefaultValue3)(const Behaviour&,
                                              const std::string&) =
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
  // wrapping the FiniteStrainBehaviourOptions::StressMeasure enum
  boost::python::enum_<FiniteStrainBehaviourOptions::StressMeasure>("FiniteStrainBehaviourOptionsStressMeasure")
    .value("CAUCHY",FiniteStrainBehaviourOptions::CAUCHY)
    .value("PK1",FiniteStrainBehaviourOptions::PK1)
    .value("PK2",FiniteStrainBehaviourOptions::PK2)
    ;
  // wrapping the FiniteStrainBehaviourOptions::TangentOperator enum
  boost::python::enum_<FiniteStrainBehaviourOptions::TangentOperator>("FiniteStrainBehaviourOptionsTangentOperator")
    .value("DSIG_DF",FiniteStrainBehaviourOptions::DSIG_DF)
    .value("DCAUCHY_DF",FiniteStrainBehaviourOptions::DSIG_DF)
    .value("DPK1_DF",FiniteStrainBehaviourOptions::DPK1_DF)
    .value("DS_DEGL",FiniteStrainBehaviourOptions::DS_DEGL)
    ;
  // wrapping the FiniteStrainBehaviourOptions class
  boost::python::class_<FiniteStrainBehaviourOptions>("FiniteStrainBehaviourOptions")
    .def_readwrite("stress_measure", &FiniteStrainBehaviourOptions::stress_measure,
		   "defines the stress measure")
    .def_readwrite("tangent_operator", &FiniteStrainBehaviourOptions::tangent_operator,
		   "defines the tangent operator")
    ;
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
      .def("getSymmetry", &Behaviour_getSymmetry,
           "return the behaviour symmetry")
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
      .add_property("params", &Behaviour_getParameters, "list of parameters")
      .add_property("parameters", &Behaviour_getParameters,
                    "list of parameters (same as the `params` property")
      .add_property("iparams", &Behaviour_getIntegerParameters,
                    "list of integer parameters")
      .add_property(
          "integer_parameters", &Behaviour_getIntegerParameters,
          "list of integer parameters (same as the `integer_params` property")
      .add_property("usparams", &Behaviour_getUnsignedShortParameters,
                    "list of unsigned short parameters")
      .add_property(
          "unsigned_short_parameters", &Behaviour_getUnsignedShortParameters,
          "list of unsigned short parameters (same as the `usparams` property)")
      .add_property("tangent_operator_blocks",
                    Behaviour_getTangentOperatorBlocks)
      .def("setParameter", setParameter1)
      .def("setIntegerParameter", setParameter2)
      .def("setUnsignedShortParameter", setParameter3)
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
  // wrapping free functions
  boost::python::def(
      "isStandardFiniteStrainBehaviour",
      mgis::behaviour::isStandardFiniteStrainBehaviour,
      "return if the given behaviour is a standard finite strain behaviour, "
      "i.e. is a finite strain behaviour using the standard finite strain "
      "kinematic (called F-Cauchy although the stress measure can be chosen "
      "when loading the behaviour)");
  boost::python::def("load", load_ptr);
  boost::python::def("load", load_ptr2);
  boost::python::def("setParameter", setParameter1);
  boost::python::def("setIntegerParameter", setParameter2);
  boost::python::def("setUnsignedShortParameter", setParameter3);
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
