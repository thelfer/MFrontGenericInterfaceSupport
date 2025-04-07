/*!
 * \file   Integrate.cxx
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
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Python/NumPySupport.hxx"

void declareIntegrate(pybind11::module_&);

static int integrateBehaviourData1(mgis::behaviour::BehaviourData& d,
                                   const mgis::behaviour::Behaviour& b) {
  auto v = mgis::behaviour::make_view(d);
  const auto s = mgis::behaviour::integrate(v, b);
  return s;
}  // end of integrateBehaviourData

static int BehaviourDataView_executeInitializeFunction(
    mgis::behaviour::BehaviourDataView& v,
    const mgis::behaviour::Behaviour& b,
    const std::string& n) {
  return mgis::behaviour::executeInitializeFunction(v, b, n);
}

static int BehaviourDataView_executeInitializeFunction2(
    mgis::behaviour::BehaviourDataView& v,
    const mgis::behaviour::Behaviour& b,
    const std::string& n,
    pybind11::object i) {
  const auto inputs = mgis::python::mgis_convert_to_span(i);
  return mgis::behaviour::executeInitializeFunction(v, b, n, inputs);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction(
    mgis::behaviour::MaterialDataManager& m, const std::string& n) {
  return mgis::behaviour::executeInitializeFunction(m, n);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction1(
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n,
    const mgis::size_type b,
    const mgis::size_type e) {
  return mgis::behaviour::executeInitializeFunction(m, n, b, e);
}

static mgis::behaviour::MultiThreadedBehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction2(
    mgis::ThreadPool& t,
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n) {
  return mgis::behaviour::executeInitializeFunction(t, m, n);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction3(
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n,
    pybind11::object i) {
  const auto inputs = mgis::python::mgis_convert_to_span(i);
  return mgis::behaviour::executeInitializeFunction(m, n, inputs);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction4(
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n,
    pybind11::object i,
    const mgis::size_type b,
    const mgis::size_type e) {
  const auto inputs = mgis::python::mgis_convert_to_span(i);
  return mgis::behaviour::executeInitializeFunction(m, n, inputs, b, e);
}

static mgis::behaviour::MultiThreadedBehaviourIntegrationResult
MaterialDataManager_executeInitializeFunction5(
    mgis::ThreadPool& t,
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n,
    pybind11::object i) {
  const auto inputs = mgis::python::mgis_convert_to_span(i);
  return mgis::behaviour::executeInitializeFunction(t, m, n, inputs);
}

static int BehaviourDataView_executePostProcessing(
    pybind11::object o,
    mgis::behaviour::BehaviourDataView& v,
    const mgis::behaviour::Behaviour& b,
    const std::string& n) {
  const auto output = mgis::python::mgis_convert_to_span(o);
  return mgis::behaviour::executePostProcessing(output, v, b, n);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executePostProcessing(
    pybind11::object o,
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n) {
  const auto output = mgis::python::mgis_convert_to_span(o);
  return mgis::behaviour::executePostProcessing(output, m, n);
}

static mgis::behaviour::BehaviourIntegrationResult
MaterialDataManager_executePostProcessing2(
    pybind11::object o,
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n,
    const mgis::size_type b,
    const mgis::size_type e) {
  const auto output = mgis::python::mgis_convert_to_span(o);
  return mgis::behaviour::executePostProcessing(output, m, n, b, e);
}

static mgis::behaviour::MultiThreadedBehaviourIntegrationResult
MaterialDataManager_executePostProcessing3(
    pybind11::object o,
    mgis::ThreadPool& t,
    mgis::behaviour::MaterialDataManager& m,
    const std::string& n) {
  const auto output = mgis::python::mgis_convert_to_span(o);
  return mgis::behaviour::executePostProcessing(output, t, m, n);
}

void declareIntegrate(pybind11::module_& m) {
  using namespace mgis::behaviour;

  pybind11::enum_<IntegrationType>(m, "IntegrationType")
      .value("PREDICTION_TANGENT_OPERATOR",
             IntegrationType::PREDICTION_TANGENT_OPERATOR)
      .value("PredictionWithTangentOperator",
             IntegrationType::PREDICTION_TANGENT_OPERATOR)
      .value("PREDICTION_SECANT_OPERATOR",
             IntegrationType::PREDICTION_SECANT_OPERATOR)
      .value("PredictionWithSecantOperator",
             IntegrationType::PREDICTION_SECANT_OPERATOR)
      .value("PREDICTION_ELASTIC_OPERATOR",
             IntegrationType::PREDICTION_ELASTIC_OPERATOR)
      .value("PredictionWithElasticOperator",
             IntegrationType::PREDICTION_ELASTIC_OPERATOR)
      .value("INTEGRATION_NO_TANGENT_OPERATOR",
             IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR)
      .value("IntegrationWithoutTangentOperator",
             IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR)
      .value("INTEGRATION_ELASTIC_OPERATOR",
             IntegrationType::INTEGRATION_ELASTIC_OPERATOR)
      .value("IntegrationWithElasticOperator",
             IntegrationType::INTEGRATION_ELASTIC_OPERATOR)
      .value("INTEGRATION_SECANT_OPERATOR",
             IntegrationType::INTEGRATION_SECANT_OPERATOR)
      .value("IntegrationWithSecantOperator",
             IntegrationType::INTEGRATION_SECANT_OPERATOR)
      .value("INTEGRATION_TANGENT_OPERATOR",
             IntegrationType::INTEGRATION_TANGENT_OPERATOR)
      .value("IntegrationWithTangentOperator",
             IntegrationType::INTEGRATION_TANGENT_OPERATOR)
      .value("INTEGRATION_CONSISTENT_TANGENT_OPERATOR",
             IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR)
      .value("IntegrationWithConsistentTangentOperator",
             IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);

  pybind11::class_<BehaviourIntegrationOptions>(m,
                                                "BehaviourIntegrationOptions")
      .def_readonly("integration_type",
                    &BehaviourIntegrationOptions::integration_type)
      .def_readonly("compute_speed_of_sound",
                    &BehaviourIntegrationOptions::compute_speed_of_sound);

  pybind11::class_<BehaviourIntegrationResult>(m, "BehaviourIntegrationResult")
      .def_readonly("exit_status", &BehaviourIntegrationResult::exit_status,
                    "The returned value has the following meaning:\n"
                    "- -1: integration failed for at least one Gauss point\n"
                    "-  0: all integrations succeeded but results are\n "
                    "      unreliable for at least one Gauss point\n"
                    "-  1: integration succeeded and results are reliable.")
      .def_readonly("time_step_increase_factor",
                    &BehaviourIntegrationResult::time_step_increase_factor)
      .def_readonly("n", &BehaviourIntegrationResult::n,
                    "number of the integration point that failed\n"
                    "or number of the last integration point \n"
                    "that reported unreliable results")
      .def_readonly("error_message",
                    &BehaviourIntegrationResult::error_message);

  pybind11::class_<MultiThreadedBehaviourIntegrationResult>(
      m, "MultiThreadedBehaviourIntegrationResult")
      .def_readonly("exit_status",
                    &MultiThreadedBehaviourIntegrationResult::exit_status,
                    "The returned value has the following meaning:\n"
                    "- -1: integration failed for at least one Gauss point\n"
                    "-  0: all integrations succeeded but results are\n "
                    "      unreliable for at least one Gauss point\n"
                    "-  1: integration succeeded and results are reliable.")
      .def_readonly("results",
                    &MultiThreadedBehaviourIntegrationResult::results);

  m.def("executeInitializeFunction",
        BehaviourDataView_executeInitializeFunction);
  m.def("executeInitializeFunction",
        BehaviourDataView_executeInitializeFunction2);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction1);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction2);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction3);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction4);
  m.def("executeInitializeFunction",
        MaterialDataManager_executeInitializeFunction5);

  int (*integrate_ptr1)(BehaviourDataView&, const Behaviour&) = integrate;
  int (*integrate_ptr2)(MaterialDataManager&, const IntegrationType,
                        const mgis::real, const mgis::size_type,
                        const mgis::size_type) = integrate;
  int (*integrate_ptr3)(mgis::ThreadPool&, MaterialDataManager&,
                        const IntegrationType, const mgis::real) = integrate;
  BehaviourIntegrationResult (*integrate_ptr4)(
      MaterialDataManager&, const BehaviourIntegrationOptions&,
      const mgis::real, const mgis::size_type, const mgis::size_type) =
      integrate;
  MultiThreadedBehaviourIntegrationResult (*integrate_ptr5)(
      mgis::ThreadPool&, MaterialDataManager&,
      const BehaviourIntegrationOptions&, const mgis::real) = integrate;

  m.def("integrate", &integrateBehaviourData1);
  m.def("integrate", integrate_ptr1);
  m.def("integrate", integrate_ptr2);
  m.def("integrate", integrate_ptr3);
  m.def("integrate", integrate_ptr4);
  m.def("integrate", integrate_ptr5);
  m.def("executePostProcessing", BehaviourDataView_executePostProcessing);
  m.def("executePostProcessing", MaterialDataManager_executePostProcessing);
  m.def("executePostProcessing", MaterialDataManager_executePostProcessing2);
  m.def("executePostProcessing", MaterialDataManager_executePostProcessing3);

}  // end of declareIntegrate
