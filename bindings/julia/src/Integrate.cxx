/*!
 * \file   bindings/julia/src/Integrate.cxx
 * \brief
 * \author Thomas Helfer
 * \date   19/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareIntegrate(jlcxx::Module& m);

static int integrateBehaviourData1(mgis::behaviour::BehaviourData& d,
                                   const mgis::behaviour::Behaviour& b) {
  auto v = mgis::behaviour::make_view(d);
  const auto s = mgis::behaviour::integrate(v, b);
  d.rdt = v.rdt;
  return s;
}  // end of integrateBehaviourData

void declareIntegrate(jlcxx::Module& m) {
  using mgis::behaviour::IntegrationType;
  m.add_bits<IntegrationType>("IntegrationType");
  m.set_const("PREDICTION_TANGENT_OPERATOR",
              IntegrationType::PREDICTION_TANGENT_OPERATOR);
  m.set_const("PredictionWithTangentOperator",
              IntegrationType::PREDICTION_TANGENT_OPERATOR);
  m.set_const("PREDICTION_SECANT_OPERATOR",
              IntegrationType::PREDICTION_SECANT_OPERATOR);
  m.set_const("PredictionWithSecantOperator",
              IntegrationType::PREDICTION_SECANT_OPERATOR);
  m.set_const("PREDICTION_ELASTIC_OPERATOR",
              IntegrationType::PREDICTION_ELASTIC_OPERATOR);
  m.set_const("PredictionWithElasticOperator",
              IntegrationType::PREDICTION_ELASTIC_OPERATOR);
  m.set_const("INTEGRATION_NO_TANGENT_OPERATOR",
              IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR);
  m.set_const("IntegrationWithoutTangentOperator",
              IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR);
  m.set_const("INTEGRATION_ELASTIC_OPERATOR",
              IntegrationType::INTEGRATION_ELASTIC_OPERATOR);
  m.set_const("IntegrationWithElasticOperator",
              IntegrationType::INTEGRATION_ELASTIC_OPERATOR);
  m.set_const("INTEGRATION_SECANT_OPERATOR",
              IntegrationType::INTEGRATION_SECANT_OPERATOR);
  m.set_const("IntegrationWithSecantOperator",
              IntegrationType::INTEGRATION_SECANT_OPERATOR);
  m.set_const("INTEGRATION_TANGENT_OPERATOR",
              IntegrationType::INTEGRATION_TANGENT_OPERATOR);
  m.set_const("IntegrationWithTangentOperator",
              IntegrationType::INTEGRATION_TANGENT_OPERATOR);
  m.set_const("INTEGRATION_CONSISTENT_TANGENT_OPERATOR",
              IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);
  m.set_const("IntegrationWithConsistentTangentOperator",
              IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);

  int (*integrate_ptr1)(mgis::behaviour::BehaviourDataView&,
                        const mgis::behaviour::Behaviour&) =
      mgis::behaviour::integrate;
  //   int (*integrate_ptr2)(mgis::behaviour::MaterialDataManager&,
  //                         const IntegrationType, const mgis::real,
  //                         const mgis::size_type, const mgis::size_type) =
  //       mgis::behaviour::integrate;
  //   int (*integrate_ptr3)(
  //       mgis::ThreadPool&, mgis::behaviour::MaterialDataManager&,
  //       const IntegrationType, const mgis::real) =
  //       mgis::behaviour::integrate;

  m.method("integrate", &integrateBehaviourData1);
  m.method("integrate", integrate_ptr1);
  //   boost::python::def("integrate", integrate_ptr2);
  //   boost::python::def("integrate", integrate_ptr3);
}  // end of declareIntegrate
