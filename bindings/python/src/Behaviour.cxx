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
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

// forward declaration
void declareBehaviour();

static void rotate_gradients_in_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &g,
    boost::python::object &r) {
  mgis::behaviour::rotateGradients(mgis::python::mgis_convert_to_span(g), b,
                                   mgis::python::mgis_convert_to_span(r));
}  // end of rotate_gradients_in_place_member

static void rotate_gradients_in_place(boost::python::object &g,
                                      const mgis::behaviour::Behaviour &b,
                                      boost::python::object &r) {
  mgis::behaviour::rotateGradients(mgis::python::mgis_convert_to_span(g), b,
                                   mgis::python::mgis_convert_to_span(r));
}  // end of rotate_gradients_in_place

static void rotate_gradients_out_of_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &mg,
    boost::python::object &gg,
    boost::python::object &r) {
  mgis::behaviour::rotateGradients(mgis::python::mgis_convert_to_span(mg), b,
                                   mgis::python::mgis_convert_to_span(gg),
                                   mgis::python::mgis_convert_to_span(r));
}  // end of rotate_gradients_out_of_place_member

static void rotate_gradients_out_of_place(boost::python::object &mg,
                                          const mgis::behaviour::Behaviour &b,
                                          boost::python::object &gg,
                                          boost::python::object &r) {
  mgis::behaviour::rotateGradients(mgis::python::mgis_convert_to_span(mg), b,
                                   mgis::python::mgis_convert_to_span(gg),
                                   mgis::python::mgis_convert_to_span(r));
}  // end of rotate_gradients_out_of_place

static void rotate_thermodynamic_forces_in_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &g,
    boost::python::object &r) {
  mgis::behaviour::rotateThermodynamicForces(
      mgis::python::mgis_convert_to_span(g), b,
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_thermodynamic_forces_in_place_member

static void rotate_thermodynamic_forces_in_place(
    boost::python::object &g,
    const mgis::behaviour::Behaviour &b,
    boost::python::object &r) {
  mgis::behaviour::rotateThermodynamicForces(
      mgis::python::mgis_convert_to_span(g), b,
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_thermodynamic_forces_in_place

static void rotate_thermodynamic_forces_out_of_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &mg,
    boost::python::object &gg,
    boost::python::object &r) {
  mgis::behaviour::rotateThermodynamicForces(
      mgis::python::mgis_convert_to_span(mg), b,
      mgis::python::mgis_convert_to_span(gg),
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_thermodynamic_forces_out_of_place_member

static void rotate_thermodynamic_forces_out_of_place(
    boost::python::object &mg,
    const mgis::behaviour::Behaviour &b,
    boost::python::object &gg,
    boost::python::object &r) {
  mgis::behaviour::rotateThermodynamicForces(
      mgis::python::mgis_convert_to_span(mg), b,
      mgis::python::mgis_convert_to_span(gg),
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_thermodynamic_forces_out_of_place

static void rotate_tangent_operator_blocks_in_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &g,
    boost::python::object &r) {
  mgis::behaviour::rotateTangentOperatorBlocks(
      mgis::python::mgis_convert_to_span(g), b,
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_tangent_operator_blocks_in_place_member

static void rotate_tangent_operator_blocks_in_place(
    boost::python::object &g,
    const mgis::behaviour::Behaviour &b,
    boost::python::object &r) {
  mgis::behaviour::rotateTangentOperatorBlocks(
      mgis::python::mgis_convert_to_span(g), b,
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_tangent_operator_blocks_in_place

static void rotate_tangent_operator_blocks_out_of_place_member(
    const mgis::behaviour::Behaviour &b,
    boost::python::object &mg,
    boost::python::object &gg,
    boost::python::object &r) {
  mgis::behaviour::rotateTangentOperatorBlocks(
      mgis::python::mgis_convert_to_span(mg), b,
      mgis::python::mgis_convert_to_span(gg),
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_tangent_operator_blocks_out_of_place_member

static void rotate_tangent_operator_blocks_out_of_place(
    boost::python::object &mg,
    const mgis::behaviour::Behaviour &b,
    boost::python::object &gg,
    boost::python::object &r) {
  mgis::behaviour::rotateTangentOperatorBlocks(
      mgis::python::mgis_convert_to_span(mg), b,
      mgis::python::mgis_convert_to_span(gg),
      mgis::python::mgis_convert_to_span(r));
}  // end of rotate_tangent_operator_blocks_out_of_place

static boost::python::list Behaviour_getInitializeFunctionsNames(
    const mgis::behaviour::Behaviour &b) {
  auto names = std::vector<std::string>{};
  for (const auto &ifct : b.initialize_functions) {
    names.push_back(ifct.first);
  }
  return mgis::python::convert_vector_to_list(names);
}  // edn of Behaviour_getInitializeFunctionsNames

static std::vector<mgis::behaviour::Variable>
Behaviour_getInitializeFunctionInputs(const mgis::behaviour::Behaviour &b,
                                      const std::string &n) {
  const auto p = b.initialize_functions.find(n);
  if (p == b.initialize_functions.end()) {
    mgis::raise(
        "getInitializeFunctionInputs: "
        "no initialize function named '" +
        n + "'");
  }
  return p->second.inputs;
}  // end of Behaviour_getInitializeFunctionInputs

static boost::python::list Behaviour_getPostProcessingsNames(
    const mgis::behaviour::Behaviour &b) {
  auto names = std::vector<std::string>{};
  for (const auto &p : b.postprocessings) {
    names.push_back(p.first);
  }
  return mgis::python::convert_vector_to_list(names);
}  // edn of Behaviour_getPostProcessingsNames

static std::vector<mgis::behaviour::Variable>
Behaviour_getPostProcessingOutputs(const mgis::behaviour::Behaviour &b,
                                   const std::string &n) {
  const auto p = b.postprocessings.find(n);
  if (p == b.postprocessings.end()) {
    mgis::raise(
        "getPostProcessingOutputs: "
        "no initialize function named '" +
        n + "'");
  }
  return p->second.outputs;
}  // end of Behaviour_getPostProcessingOutputs

void declareBehaviour() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::BehaviourDescription;
  using mgis::behaviour::FiniteStrainBehaviourOptions;
  using mgis::behaviour::Hypothesis;
  Behaviour (*load_ptr)(const std::string &, const std::string &,
                        const Hypothesis) = &mgis::behaviour::load;
  Behaviour (*load_ptr2)(const FiniteStrainBehaviourOptions &,
                         const std::string &, const std::string &,
                         const Hypothesis) = &mgis::behaviour::load;
  void (*setParameter1)(const Behaviour &, const std::string &, const double) =
      &mgis::behaviour::setParameter;
  void (*setParameter2)(const Behaviour &, const std::string &, const int) =
      &mgis::behaviour::setParameter;
  void (*setParameter3)(const Behaviour &, const std::string &,
                        const unsigned short) = &mgis::behaviour::setParameter;
  // wrapping the FiniteStrainBehaviourOptions::StressMeasure enum
  boost::python::enum_<FiniteStrainBehaviourOptions::StressMeasure>(
      "FiniteStrainBehaviourOptionsStressMeasure")
      .value("CAUCHY", FiniteStrainBehaviourOptions::CAUCHY)
      .value("PK1", FiniteStrainBehaviourOptions::PK1)
      .value("PK2", FiniteStrainBehaviourOptions::PK2);
  // wrapping the FiniteStrainBehaviourOptions::TangentOperator enum
  boost::python::enum_<FiniteStrainBehaviourOptions::TangentOperator>(
      "FiniteStrainBehaviourOptionsTangentOperator")
      .value("DSIG_DF", FiniteStrainBehaviourOptions::DSIG_DF)
      .value("DCAUCHY_DF", FiniteStrainBehaviourOptions::DSIG_DF)
      .value("DPK1_DF", FiniteStrainBehaviourOptions::DPK1_DF)
      .value("DS_DEGL", FiniteStrainBehaviourOptions::DS_DEGL)
      .value("DTAU_DDF", FiniteStrainBehaviourOptions::DTAU_DDF);
  // wrapping the FiniteStrainBehaviourOptions class
  boost::python::class_<FiniteStrainBehaviourOptions>(
      "FiniteStrainBehaviourOptions")
      .def_readwrite("stress_measure",
                     &FiniteStrainBehaviourOptions::stress_measure,
                     "defines the stress measure")
      .def_readwrite("tangent_operator",
                     &FiniteStrainBehaviourOptions::tangent_operator,
                     "defines the tangent operator");
  // wrapping the Behaviour class
  boost::python::class_<Behaviour, boost::python::bases<BehaviourDescription>>(
      "Behaviour")
      // boost::python::class_<Behaviour>("Behaviour")
      .def("getInitializeFunctionsNames", Behaviour_getInitializeFunctionsNames)
      .def("getInitializeFunctionInputs", Behaviour_getInitializeFunctionInputs)
      .def("getPostProcessingsNames", Behaviour_getPostProcessingsNames)
      .def("getPostProcessingOutputs", Behaviour_getPostProcessingOutputs)
      .def("setParameter", setParameter1)
      .def("setIntegerParameter", setParameter2)
      .def("setUnsignedShortParameter", setParameter3)
      .def("rotateGradients", rotate_gradients_in_place_member)
      .def("rotateGradients", rotate_gradients_out_of_place_member)
      .def("rotateThermodynamicForces",
           rotate_thermodynamic_forces_in_place_member)
      .def("rotateThermodynamicForces",
           rotate_thermodynamic_forces_out_of_place_member)
      .def("rotateTangentOperatorBlocks",
           rotate_tangent_operator_blocks_in_place_member)
      .def("rotateTangentOperatorBlocks",
           rotate_tangent_operator_blocks_out_of_place_member);
  // wrapping free functions
  boost::python::def("rotateGradients", rotate_gradients_in_place);
  boost::python::def("rotateGradients", rotate_gradients_out_of_place);
  boost::python::def("rotateThermodynamicForces",
                     rotate_thermodynamic_forces_in_place);
  boost::python::def("rotateThermodynamicForces",
                     rotate_thermodynamic_forces_out_of_place);
  boost::python::def("rotateTangentOperatorBlocks",
                     rotate_tangent_operator_blocks_in_place);
  boost::python::def("rotateTangentOperatorBlocks",
                     rotate_tangent_operator_blocks_out_of_place);

  boost::python::def("load", load_ptr);
  boost::python::def("load", load_ptr2);
  boost::python::def("setParameter", setParameter1);
  boost::python::def("setIntegerParameter", setParameter2);
  boost::python::def("setUnsignedShortParameter", setParameter3);

}  // end of declareBehaviour
