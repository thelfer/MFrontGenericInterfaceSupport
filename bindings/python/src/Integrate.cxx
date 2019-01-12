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

#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

void declareIntegrate();

static int integrateBehaviourData1(mgis::behaviour::BehaviourData& d,
                                   const mgis::behaviour::Behaviour& b) {
  auto v = mgis::behaviour::make_view(d);
  const auto s = mgis::behaviour::integrate(v, b);
  d.rdt = v.rdt;
  return s;
}  // end of integrateBehaviourData

void declareIntegrate() {
  boost::python::enum_<mgis::behaviour::IntegrationType>("IntegrationType")
      .value("PREDICTION_TANGENT_OPERATOR",
             mgis::behaviour::IntegrationType::PREDICTION_TANGENT_OPERATOR)
      .value("PredictionWithTangentOperator",
             mgis::behaviour::IntegrationType::PREDICTION_TANGENT_OPERATOR)
      .value("PREDICTION_SECANT_OPERATOR",
             mgis::behaviour::IntegrationType::PREDICTION_SECANT_OPERATOR)
      .value("PredictionWithSecantOperator",
             mgis::behaviour::IntegrationType::PREDICTION_SECANT_OPERATOR)
      .value("PREDICTION_ELASTIC_OPERATOR",
             mgis::behaviour::IntegrationType::PREDICTION_ELASTIC_OPERATOR)
      .value("PredictionWithElasticOperator",
             mgis::behaviour::IntegrationType::PREDICTION_ELASTIC_OPERATOR)
      .value("INTEGRATION_NO_TANGENT_OPERATOR",
             mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR)
      .value("IntegrationWithoutTangentOperator",
             mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR)
      .value("INTEGRATION_ELASTIC_OPERATOR",
             mgis::behaviour::IntegrationType::INTEGRATION_ELASTIC_OPERATOR)
      .value("IntegrationWithElasticOperator",
             mgis::behaviour::IntegrationType::INTEGRATION_ELASTIC_OPERATOR)
      .value("INTEGRATION_SECANT_OPERATOR",
             mgis::behaviour::IntegrationType::INTEGRATION_SECANT_OPERATOR)
      .value("IntegrationWithSecantOperator",
             mgis::behaviour::IntegrationType::INTEGRATION_SECANT_OPERATOR)
      .value("INTEGRATION_TANGENT_OPERATOR",
             mgis::behaviour::IntegrationType::INTEGRATION_TANGENT_OPERATOR)
      .value("IntegrationWithTangentOperator",
             mgis::behaviour::IntegrationType::INTEGRATION_TANGENT_OPERATOR)
      .value("INTEGRATION_CONSISTENT_TANGENT_OPERATOR",
             mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR)
      .value("IntegrationWithConsistentTangentOperator",
             mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);

  int (*integrate_ptr1)(mgis::behaviour::BehaviourDataView&,
                        const mgis::behaviour::Behaviour&) =
      mgis::behaviour::integrate;
  int (*integrate_ptr2)(mgis::behaviour::MaterialDataManager&,
                        const mgis::behaviour::IntegrationType,
                        const mgis::real, const mgis::size_type,
                        const mgis::size_type) = mgis::behaviour::integrate;
  int (*integrate_ptr3)(mgis::ThreadPool&,
                        mgis::behaviour::MaterialDataManager&,
                        const mgis::behaviour::IntegrationType,
                        const mgis::real) = mgis::behaviour::integrate;

  boost::python::def("integrate", &integrateBehaviourData1);
  boost::python::def("integrate", integrate_ptr1);
  boost::python::def("integrate", integrate_ptr2);
  boost::python::def("integrate", integrate_ptr3);
}  // end of declareIntegrate
