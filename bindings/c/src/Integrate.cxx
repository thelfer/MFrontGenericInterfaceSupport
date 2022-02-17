/*!
 * \file   Integrate.cxx
 * \brief
 * \author th202608
 * \date   02/08/2018
 * \copyright Copyright (C) 2006-2018 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/BehaviourDataView.h"
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/Integrate.h"
#include "MGIS/Behaviour/Integrate.hxx"

extern "C" {

static mgis::behaviour::IntegrationType convertIntegrationType(
    const mgis_bv_IntegrationType i) {
  switch (i) {
    case MGIS_BV_PREDICTION_TANGENT_OPERATOR:
      return mgis::behaviour::IntegrationType::PREDICTION_TANGENT_OPERATOR;
    case MGIS_BV_PREDICTION_SECANT_OPERATOR:
      return mgis::behaviour::IntegrationType::PREDICTION_SECANT_OPERATOR;
    case MGIS_BV_PREDICTION_ELASTIC_OPERATOR:
      return mgis::behaviour::IntegrationType::PREDICTION_ELASTIC_OPERATOR;
    case MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR:
      return mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
    case MGIS_BV_INTEGRATION_ELASTIC_OPERATOR:
      return mgis::behaviour::IntegrationType::INTEGRATION_ELASTIC_OPERATOR;
    case MGIS_BV_INTEGRATION_SECANT_OPERATOR:
      return mgis::behaviour::IntegrationType::INTEGRATION_SECANT_OPERATOR;
    case MGIS_BV_INTEGRATION_TANGENT_OPERATOR:
      return mgis::behaviour::IntegrationType::INTEGRATION_TANGENT_OPERATOR;
    case MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR:
      return mgis::behaviour::IntegrationType::
          INTEGRATION_CONSISTENT_TANGENT_OPERATOR;
    default:
      mgis::raise("convertIntegrationType: invalid integration type");
  }
  return mgis::behaviour::IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
}  // end of convertIntegrationType

mgis_status mgis_bv_integrate(int* const r,
                              mgis_bv_BehaviourDataView* const d,
                              const mgis_bv_Behaviour* const b) {
  *r = mgis::behaviour::integrate(*d, *b);
  if ((*r != 1) && (*r != 0)) {
    return mgis_report_failure("behaviour integration failed");
  }
  return mgis_report_success();
}  // end of mgis_bv_integrate

mgis_status mgis_bv_integrate_2(int* const r,
                                mgis_bv_BehaviourData* const d,
                                const mgis_bv_Behaviour* const b) {
  auto v = mgis::behaviour::make_view(*d);
  auto s = mgis_bv_integrate(r, &v, b);
  return s;
}  // end of mgis_bv_integrate2

mgis_status mgis_bv_integrate_material_data_manager(
    int* const r,
    mgis_ThreadPool* const p,
    mgis_bv_MaterialDataManager* const m,
    const mgis_bv_IntegrationType i,
    const mgis_real dt) {
  *r = -1;
  try {
    *r = mgis::behaviour::integrate(*p, *m, convertIntegrationType(i), dt);
    if ((*r != 1) && (*r != 0)) {
      return mgis_report_failure("behaviour integration failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_integrate_material_data_manager

mgis_status mgis_bv_integrate_material_data_manager_part(
    int* const r,
    mgis_bv_MaterialDataManager* const m,
    const mgis_bv_IntegrationType i,
    const mgis_real dt,
    const mgis_size_type b,
    const mgis_size_type e) {
  *r = -1;
  try {
    *r = mgis::behaviour::integrate(*m, convertIntegrationType(i), dt, b, e);
    if ((*r != 1) && (*r != 0)) {
      return mgis_report_failure("behaviour integration failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_integrate_material_data_manager_part

}  // end of extern "C"
