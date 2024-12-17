/*!
 * \file   Behaviour.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <limits>
#include <cstring>
#include <utility>
#include <stdexcept>
#include "MGIS/Behaviour/Behaviour.h"

extern "C" {

mgis_status mgis_bv_create_finite_strain_behaviour_options(
    mgis_bv_FiniteStrainBehaviourOptions** o) {
  *o = nullptr;
  try {
    *o = new mgis::behaviour::FiniteStrainBehaviourOptions();
    if (*o == nullptr) {
      return mgis_report_failure(
          "mgis_bv_create_finite_strain_behaviour_options: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_finite_strain_behaviour_options

mgis_status mgis_bv_finite_strain_behaviour_options_set_stress_measure(
    mgis_bv_FiniteStrainBehaviourOptions* const o,
    const mgis_bv_FiniteStrainBehavourStressMeasure s) {
  if (o == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  auto sm = mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY;
  if (s == MGIS_BV_CAUCHY) {
  } else if (s == MGIS_BV_PK2) {
    sm = mgis::behaviour::FiniteStrainBehaviourOptions::PK2;
  } else if (s == MGIS_BV_PK1) {
    sm = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
  } else {
    return mgis_report_failure("invalid stress measure");
  }
  o->stress_measure = sm;
  return mgis_report_success();
}  // end of mgis_bv_finite_strain_behaviour_options_set_stress_measure

mgis_status
mgis_bv_finite_strain_behaviour_options_set_stress_measure_by_string(
    mgis_bv_FiniteStrainBehaviourOptions* const o, const char* const s) {
  if ((o == nullptr) || (s == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  auto sm = mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY;
  if ((strcmp(s, "CauchyStress") == 0) || (strcmp(s, "CAUCHY") == 0)) {
  } else if ((strcmp(s, "SecondPiolaKirchhoffStress") == 0) ||
             (strcmp(s, "PK2") == 0) || (strcmp(s, "S") == 0)) {
    sm = mgis::behaviour::FiniteStrainBehaviourOptions::PK2;
  } else if ((strcmp(s, "FirstPiolaKirchhoffStress") == 0) ||
             (strcmp(s, "PK1") == 0)) {
    sm = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
  } else {
    return mgis_report_failure("invalid stress measure");
  }
  o->stress_measure = sm;
  return mgis_report_success();
}  // end of
   // mgis_bv_finite_strain_behaviour_options_set_stress_measure_by_string

mgis_status mgis_bv_finite_strain_behaviour_options_set_tangent_operator(
    mgis_bv_FiniteStrainBehaviourOptions* const o,
    const mgis_bv_FiniteStrainBehavourTangentOperator t) {
  if (o == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  auto tm = mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF;
  if (t == MGIS_BV_DSIG_DF) {
  } else if (t == MGIS_BV_DS_DEGL) {
    tm = mgis::behaviour::FiniteStrainBehaviourOptions::DS_DEGL;
  } else if (t == MGIS_BV_DPK1_DF) {
    tm = mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF;
  } else {
    return mgis_report_failure("invalid tangent operator");
  }
  o->tangent_operator = tm;
  return mgis_report_success();
}  // end of
   // mgis_bv_finite_strain_behaviour_options_set_tangent_operator

mgis_status
mgis_bv_finite_strain_behaviour_options_set_tangent_operator_by_string(
    mgis_bv_FiniteStrainBehaviourOptions* const o, const char* const t) {
  if ((o == nullptr) || (t == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  auto tm = mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF;
  if ((strcmp(t, "DSIG_DF") == 0) || (strcmp(t, "DCAUCHY_DF") == 0)) {
  } else if (strcmp(t, "DS_DEGL") == 0) {
    tm = mgis::behaviour::FiniteStrainBehaviourOptions::DS_DEGL;
  } else if (strcmp(t, "DPK1_DF") == 0) {
    tm = mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF;
  } else {
    return mgis_report_failure("invalid tangent operator");
  }
  o->tangent_operator = tm;
  return mgis_report_success();
}  // end of
   // mgis_bv_finite_strain_behaviour_options_set_tangent_operator_by_string

mgis_status mgis_bv_free_finite_strain_behaviour_options(
    mgis_bv_FiniteStrainBehaviourOptions** o) {
  try {
    delete *o;
    *o = nullptr;
  } catch (...) {
    *o = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_free_finite_strain_behaviour_options

mgis_status mgis_bv_is_standard_finite_strain_behaviour(int* r,
                                                        const char* const l,
                                                        const char* const b) {
  *r = static_cast<int>(false);
  try {
    *r = static_cast<int>(
        mgis::behaviour::isStandardFiniteStrainBehaviour(l, b));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_is_standard_finite_strain_behaviour

mgis_status mgis_bv_load_behaviour(mgis_bv_Behaviour** ptr,
                                   const char* const l,
                                   const char* const b,
                                   const char* const h) {
  *ptr = nullptr;
  try {
    const auto bv = mgis::behaviour::load(l, b, mgis::behaviour::fromString(h));
    *ptr = new mgis::behaviour::Behaviour(std::move(bv));
    if (*ptr == nullptr) {
      return mgis_report_failure(
          "mgis_bv_load_behaviour: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_load_behaviour

mgis_status mgis_bv_load_finite_strain_behaviour(
    mgis_bv_Behaviour** ptr,
    const mgis_bv_FiniteStrainBehaviourOptions* const o,
    const char* const l,
    const char* const b,
    const char* const h) {
  *ptr = nullptr;
  try {
    const auto bv =
        mgis::behaviour::load(*o, l, b, mgis::behaviour::fromString(h));
    *ptr = new mgis::behaviour::Behaviour(std::move(bv));
    if (*ptr == nullptr) {
      return mgis_report_failure(
          "mgis_bv_load_behaviour: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_load_finite_strain_behaviour

mgis_status mgis_bv_behaviour_rotate_gradients_in_place(
    mgis_real* const g,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r) {
  b->rotate_gradients_ptr(g, g, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_gradients_in_place

mgis_status mgis_bv_behaviour_rotate_gradients_out_of_place(
    mgis_real* const mg,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gg,
    const mgis_real* const r) {
  b->rotate_gradients_ptr(mg, gg, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_gradients_out_of_place

mgis_status mgis_bv_behaviour_rotate_array_of_gradients_in_place(
    mgis_real* const g,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_gradients_ptr(g, g, r, s);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_array_of_gradients_in_place

mgis_status mgis_bv_behaviour_rotate_array_of_gradients_out_of_place(
    mgis_real* const mg,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gg,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_gradients_ptr(mg, gg, r, s);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_array_of_gradients_out_of_place

mgis_status mgis_bv_behaviour_rotate_thermodynamic_forces_in_place(
    mgis_real* const th,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r) {
  b->rotate_thermodynamic_forces_ptr(th, th, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_thermodynamic_forces_in_place

mgis_status mgis_bv_behaviour_rotate_thermodynamic_forces_out_of_place(
    mgis_real* const mth,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gth,
    const mgis_real* const r) {
  b->rotate_thermodynamic_forces_ptr(mth, gth, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_thermodynamic_forces_out_of_place

mgis_status mgis_bv_behaviour_rotate_array_of_thermodynamic_forces_in_place(
    mgis_real* const th,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_thermodynamic_forces_ptr(th, th, r, s);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_array_of_thermodynamic_forces_in_place

mgis_status mgis_bv_behaviour_rotate_array_of_thermodynamic_forces_out_of_place(
    mgis_real* const mth,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gth,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_thermodynamic_forces_ptr(mth, gth, r, s);
  return mgis_report_success();
}  // end of
   // mgis_bv_behaviour_rotate_array_of_thermodynamic_forces_out_of_place

mgis_status mgis_bv_behaviour_rotate_tangent_operator_blocks_in_place(
    mgis_real* const to,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r) {
  b->rotate_tangent_operator_blocks_ptr(to, to, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_tangent_operator_blocks_in_place

mgis_status mgis_bv_behaviour_rotate_tangent_operator_blocks_out_of_place(
    mgis_real* const mto,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gto,
    const mgis_real* const r) {
  b->rotate_tangent_operator_blocks_ptr(mto, gto, r);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_tangent_operator_blocks_out_of_place

mgis_status mgis_bv_behaviour_rotate_array_of_tangent_operator_blocks_in_place(
    mgis_real* const to,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_tangent_operator_blocks_ptr(to, to, r, s);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_rotate_array_of_tangent_operator_blocks_in_place

mgis_status
mgis_bv_behaviour_rotate_array_of_tangent_operator_blocks_out_of_place(
    mgis_real* const mto,
    const mgis_bv_Behaviour* const b,
    const mgis_real* const gto,
    const mgis_real* const r,
    const mgis_size_type s) {
  b->rotate_array_of_tangent_operator_blocks_ptr(mto, gto, r, s);
  return mgis_report_success();
}  // end of
   // mgis_bv_behaviour_rotate_array_of_tangent_operator_blocks_out_of_place

mgis_status mgis_bv_behaviour_get_tangent_operator_array_size(
    mgis::size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = mgis::size_type{};
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = mgis::behaviour::getTangentOperatorArraySize(*b);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_tangent_operator_array_size

mgis_status mgis_bv_behaviour_get_library(const char** l,
                                          const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *l = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *l = b->library.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_library

mgis_status mgis_bv_behaviour_get_source(const char** s,
                                         const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->source.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_source

mgis_status mgis_bv_behaviour_get_author(const char** s,
                                         const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->author.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_author

mgis_status mgis_bv_behaviour_get_date(const char** s,
                                       const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->date.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_date

mgis_status mgis_bv_behaviour_get_validator(const char** s,
                                            const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->validator.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_validator

mgis_status mgis_bv_behaviour_get_build_id(const char** s,
                                           const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->build_id.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_build_id

mgis_status mgis_bv_behaviour_get_tfel_version(
    const char** v, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *v = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *v = b->tfel_version.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_tfel_version

mgis_status mgis_bv_behaviour_get_behaviour_name(
    const char** h, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = b->behaviour.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_behaviour_name

mgis_status mgis_bv_behaviour_get_function_name(
    const char** h, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = b->function.c_str();
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_function_name

mgis_status mgis_bv_behaviour_get_hypothesis(const char** h,
                                             const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = mgis::behaviour::toString(b->hypothesis);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_hypothesis

MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_symmetry(
    mgis_bv_BehaviourSymmetry* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  switch (b->symmetry) {
    case mgis::behaviour::Behaviour::ISOTROPIC:
      *s = MGIS_BV_ISOTROPIC;
      break;
    case mgis::behaviour::Behaviour::ORTHOTROPIC:
      *s = MGIS_BV_ORTHOTROPIC;
      break;
    default:
      return mgis_report_failure("unsupported behaviour symmetry type");
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_behaviour_symmetry

MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_type(
    mgis_bv_BehaviourType* const t, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  switch (b->btype) {
    case mgis::behaviour::Behaviour::GENERALBEHAVIOUR:
      *t = MGIS_BV_GENERALBEHAVIOUR;
      break;
    case mgis::behaviour::Behaviour::STANDARDSTRAINBASEDBEHAVIOUR:
      *t = MGIS_BV_STANDARDSTRAINBASEDBEHAVIOUR;
      break;
    case mgis::behaviour::Behaviour::STANDARDFINITESTRAINBEHAVIOUR:
      *t = MGIS_BV_STANDARDFINITESTRAINBEHAVIOUR;
      break;
    case mgis::behaviour::Behaviour::COHESIVEZONEMODEL:
      *t = MGIS_BV_COHESIVEZONEMODEL;
      break;
    default:
      return mgis_report_failure("unsupported behaviour type");
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_behaviour_type

MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_kinematic(
    mgis_bv_BehaviourKinematic* const k, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *k = MGIS_BV_UNDEFINEDKINEMATIC;
    return mgis_report_failure("invalid argument");
  }
  switch (b->kinematic) {
    case mgis::behaviour::Behaviour::UNDEFINEDKINEMATIC:
      *k = MGIS_BV_UNDEFINEDKINEMATIC;
      break;
    case mgis::behaviour::Behaviour::SMALLSTRAINKINEMATIC:
      *k = MGIS_BV_SMALLSTRAINKINEMATIC;
      break;
    case mgis::behaviour::Behaviour::COHESIVEZONEKINEMATIC:
      *k = MGIS_BV_COHESIVEZONEKINEMATIC;
      break;
    case mgis::behaviour::Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY:
      *k = MGIS_BV_FINITESTRAINKINEMATIC_F_CAUCHY;
      break;
    case mgis::behaviour::Behaviour::FINITESTRAINKINEMATIC_ETO_PK1:
      *k = MGIS_BV_FINITESTRAINKINEMATIC_ETO_PK1;
      break;
    default:
      return mgis_report_failure("unsupported behaviour kinematic");
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_behaviour_kinematic

mgis_status mgis_bv_behaviour_get_gradients_size(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = mgis::behaviour::getArraySize(b->gradients, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_gradients_size

mgis_status mgis_bv_behaviour_get_thermodynamic_forces_size(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = mgis::behaviour::getArraySize(b->thermodynamic_forces, b->hypothesis);
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_thermodynamic_forces_size

mgis_status mgis_bv_behaviour_get_number_of_material_properties(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  *s = 0;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->mps.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_material_properties

mgis_status mgis_bv_behaviour_get_material_property_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& mp = b->mps.at(i);
    *n = mp.name.c_str();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_material_property_name

mgis_status mgis_bv_behaviour_get_material_property_type(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& ev = b->mps.at(i);
    switch (ev.type) {
      case mgis::behaviour::Variable::SCALAR:
        *t = MGIS_BV_SCALAR;
        break;
      case mgis::behaviour::Variable::VECTOR:
        *t = MGIS_BV_VECTOR;
        break;
      case mgis::behaviour::Variable::STENSOR:
        *t = MGIS_BV_STENSOR;
        break;
      case mgis::behaviour::Variable::TENSOR:
        *t = MGIS_BV_TENSOR;
        break;
      default:
        return mgis_report_failure("unsupported variable type");
    }
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_material_property_type

mgis_status mgis_bv_behaviour_get_number_of_internal_state_variables(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  *s = 0;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->isvs.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_internal_state_variables

mgis_status mgis_bv_behaviour_get_internal_state_variable_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = b->isvs.at(i);
    *n = iv.name.c_str();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_internal_state_variable_name

namespace c_internals {

  static mgis_status getVariableType(mgis_bv_VariableType& t,
                                     const mgis::behaviour::Variable& v) {
    switch (v.type) {
      case mgis::behaviour::Variable::SCALAR:
        t = MGIS_BV_SCALAR;
        break;
      case mgis::behaviour::Variable::VECTOR:
        t = MGIS_BV_VECTOR;
        break;
      case mgis::behaviour::Variable::STENSOR:
        t = MGIS_BV_STENSOR;
        break;
      case mgis::behaviour::Variable::TENSOR:
        t = MGIS_BV_TENSOR;
        break;
      default:
        return mgis_report_failure("unsupported variable type");
    }
    return mgis_report_success();
  }

}  // end of namespace c_internals

mgis_status mgis_bv_behaviour_get_internal_state_variable_type(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    return c_internals::getVariableType(*t, b->isvs.at(i));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_internal_state_variable_type

mgis_status mgis_bv_behaviour_get_internal_state_variable_type_by_name(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = getVariable(b->isvs, n);
    return c_internals::getVariableType(*t, iv);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_internal_state_variable_type_by_name

mgis_status mgis_bv_behaviour_get_internal_state_variable_offset(
    mgis_size_type* const o,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = b->isvs.at(i);
    *o = getVariableOffset(b->isvs, iv.name, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_internal_state_variable_offset

mgis_status mgis_bv_behaviour_get_internal_state_variable_offset_by_name(
    mgis_size_type* const o,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = getVariable(b->isvs, n);
    *o = getVariableOffset(b->isvs, iv.name, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_internal_state_variable_offset_by_name

mgis_status mgis_bv_behaviour_get_internal_state_variables_size(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = getArraySize(b->isvs, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_internal_state_variables_size

mgis_status mgis_bv_behaviour_get_number_of_external_state_variables(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  *s = 0;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->esvs.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_external_state_variables

mgis_status mgis_bv_behaviour_get_external_state_variable_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& ev = b->esvs.at(i);
    *n = ev.name.c_str();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_external_state_variable_name

mgis_status mgis_bv_behaviour_get_external_state_variable_offset(
    mgis_size_type* const o,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& esv = b->esvs.at(i);
    *o = getVariableOffset(b->esvs, esv.name, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_external_state_variable_offset

mgis_status mgis_bv_behaviour_get_external_state_variable_offset_by_name(
    mgis_size_type* const o,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& esv = getVariable(b->esvs, n);
    *o = getVariableOffset(b->esvs, esv.name, b->hypothesis);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_external_state_variable_offset_by_name

mgis_status mgis_bv_behaviour_get_external_state_variable_type(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& ev = b->esvs.at(i);
    return c_internals::getVariableType(*t, ev);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_external_state_variable_type

mgis_status mgis_bv_behaviour_get_external_state_variable_type_by_name(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& ev = getVariable(b->esvs, n);
    return c_internals::getVariableType(*t, ev);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_external_state_variable_type_by_name

mgis_status mgis_bv_behaviour_get_number_of_parameters(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->params.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_parameters

mgis_status mgis_bv_behaviour_get_parameter_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& p = b->params.at(i);
    *n = p.c_str();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_parameter_name

mgis_status mgis_bv_behaviour_get_number_of_integer_parameters(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  *s = 0;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->iparams.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_integer_parameters

mgis_status mgis_bv_behaviour_get_integer_parameter_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& p = b->iparams.at(i);
    *n = p.c_str();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_integer_parameter_name

mgis_status mgis_bv_behaviour_get_number_of_unsigned_short_parameters(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  *s = 0;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->usparams.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_behaviour_get_number_of_unsigned_short_parameters

mgis_status mgis_bv_behaviour_get_unsigned_short_parameter_name(
    const char** n, const mgis_bv_Behaviour* const b, const mgis_size_type i) {
  *n = nullptr;
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& p = b->usparams.at(i);
    *n = p.c_str();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_unsigned_short_parameter_name

mgis_status mgis_bv_behaviour_set_parameter(const mgis_bv_Behaviour* const b,
                                            const char* const n,
                                            const double v) {
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    setParameter(*b, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_set_parameter

mgis_status mgis_bv_behaviour_set_integer_parameter(
    const mgis_bv_Behaviour* const b, const char* const n, const int v) {
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    mgis::behaviour::setParameter(*b, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_set_integer_parameter

mgis_status mgis_bv_behaviour_set_unsigned_short_parameter(
    const mgis_bv_Behaviour* const b,
    const char* const n,
    const unsigned short v) {
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    mgis::behaviour::setParameter(*b, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_set_unsigned_short_parameter

mgis_status mgis_bv_behaviour_get_parameter_default_value(
    double* const v, const mgis_bv_Behaviour* const b, const char* const n) {
  *v = std::numeric_limits<double>::quiet_NaN();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getParameterDefaultValue<double>(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_parameter_default_value

mgis_status mgis_bv_behaviour_get_integer_parameter_default_value(
    int* const v, const mgis_bv_Behaviour* const b, const char* const n) {
  *v = std::numeric_limits<int>::max();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getParameterDefaultValue<int>(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_integer_parameter_default_value

mgis_status mgis_bv_behaviour_get_unsigned_short_parameter_default_value(
    unsigned short* const v,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  *v = std::numeric_limits<unsigned short>::max();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getParameterDefaultValue<unsigned short>(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_unsigned_short_parameter_default_value

mgis_status mgis_bv_behaviour_has_bounds(int* const v,
                                         const mgis_bv_Behaviour* const b,
                                         const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasBounds(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_bounds

mgis_status mgis_bv_behaviour_has_lower_bound(int* const v,
                                              const mgis_bv_Behaviour* const b,
                                              const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasLowerBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_lower_bound

mgis_status mgis_bv_behaviour_has_upper_bound(int* const v,
                                              const mgis_bv_Behaviour* const b,
                                              const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasUpperBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_upper_bound

mgis_status mgis_bv_behaviour_get_lower_bound(long double* const v,
                                              const mgis_bv_Behaviour* const b,
                                              const char* const n) {
  *v = std::numeric_limits<long double>::quiet_NaN();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getLowerBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_lower_bound

mgis_status mgis_bv_behaviour_get_upper_bound(long double* const v,
                                              const mgis_bv_Behaviour* const b,
                                              const char* const n) {
  *v = std::numeric_limits<long double>::quiet_NaN();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getUpperBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_upper_bound

mgis_status mgis_bv_behaviour_has_physical_bounds(
    int* const v, const mgis_bv_Behaviour* const b, const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasPhysicalBounds(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_physical_bounds

mgis_status mgis_bv_behaviour_has_lower_physical_bound(
    int* const v, const mgis_bv_Behaviour* const b, const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasLowerPhysicalBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_lower_physical_bound

mgis_status mgis_bv_behaviour_has_upper_physical_bound(
    int* const v, const mgis_bv_Behaviour* const b, const char* const n) {
  *v = -1;
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::hasUpperPhysicalBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_has_upper_physical_bound

mgis_status mgis_bv_behaviour_get_lower_physical_bound(
    long double* const v,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  *v = std::numeric_limits<long double>::quiet_NaN();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getLowerPhysicalBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_lower_physical_bound

mgis_status mgis_bv_behaviour_get_upper_physical_bound(
    long double* const v,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  *v = std::numeric_limits<long double>::quiet_NaN();
  if ((b == nullptr) || (n == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  try {
    *v = mgis::behaviour::getUpperPhysicalBound(*b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_behaviour_get_upper_physical_bound

mgis_status mgis_bv_behaviour_computes_stored_energy(
    int* const v, const mgis_bv_Behaviour* const b) {
  if ((b == nullptr) || (v == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  *v = static_cast<int>(b->computesStoredEnergy);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_computes_stored_energy

mgis_status mgis_bv_behaviour_computes_dissipated_energy(
    int* const v, const mgis_bv_Behaviour* const b) {
  if ((b == nullptr) || (v == nullptr)) {
    return mgis_report_failure("invalid argument");
  }
  *v = static_cast<int>(b->computesDissipatedEnergy);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_computes_dissipated_energy

mgis_status mgis_bv_free_behaviour(mgis_bv_Behaviour** b) {
  try {
    delete *b;
    *b = nullptr;
  } catch (...) {
    *b = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_free_Behaviour

}  // end of extern "C"
