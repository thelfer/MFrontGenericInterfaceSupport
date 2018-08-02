/*!
 * \file   Description.cxx
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

#include <utility>
#include <stdexcept>
#include "MGIS/Behaviour/Behaviour.h"

extern "C" {

mgis_status mgis_bv_load_behaviour(mgis_bv_Behaviour** ptr,
                                     const char* const l,
                                     const char* const b,
                                     const char* const h) {
  *ptr = nullptr;
  try {
    const auto bv = mgis::behaviour::load(l, b, mgis::behaviour::fromString(h));
    *ptr = new mgis::behaviour::Behaviour(std::move(bv));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of load

mgis_status mgis_bv_get_library(const char** l,
                                const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *l = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *l = b->library.c_str();
  return mgis_report_success();
}  // end of mgis_bv_get_library

mgis_status mgis_bv_get_source(const char** s,
                               const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *s = b->source.c_str();
  return mgis_report_success();
}  // end of mgis_bv_get_source

mgis_status mgis_bv_get_tfel_version(const char** v,
                                     const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *v = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *v = b->tfel_version.c_str();
  return mgis_report_success();
}  // end of mgis_bv_get_tfel_version

mgis_status mgis_bv_get_behaviour_name(const char** h,
                                   const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = b->behaviour.c_str();
  return mgis_report_success();
}  // end of mgis_bv_get_behaviour_name

mgis_status mgis_bv_get_function_name(const char** h,
                                   const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = b->function.c_str();
  return mgis_report_success();
}  // end of mgis_bv_get_function_name

mgis_status mgis_bv_get_hypothesis(const char** h,
                                   const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *h = nullptr;
    return mgis_report_failure("invalid argument");
  }
  *h = mgis::behaviour::toString(b->hypothesis);
  return mgis_report_success();
}  // end of mgis_bv_get_hypothesis

MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_symmetry(
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
}  // end of mgis_bv_get_behaviour_symmetry

MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_type(
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
}  // end of mgis_bv_get_behaviour_type

MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_kinematic(
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
}  // end of mgis_bv_get_behaviour_kinematic

mgis_status mgis_bv_get_number_of_material_properties(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->mps.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_number_of_material_properties

mgis_status mgis_bv_get_material_property_name(
    const char** n,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    *n = nullptr;
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& mp = b->mps.at(i);
    *n = mp.name.c_str();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_material_property_name

mgis_status mgis_bv_get_number_of_internal_state_variables(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->isvs.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_number_of_internal_state_variables

mgis_status mgis_bv_get_internal_state_variable_name(
    const char** n,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    *n = nullptr;
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = b->isvs.at(i);
    *n = iv.name.c_str();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_internal_state_variable_name

mgis_status mgis_bv_get_internal_state_variable_type(
    mgis_bv_VariableType* const t,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& iv = b->isvs.at(i);
    switch(iv.type){
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
}  // end of mgis_bv_get_internal_state_variable_type

mgis_status mgis_bv_get_number_of_external_state_variables(
    mgis_size_type* const s, const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *s = 0;
    return mgis_report_failure("invalid argument");
  }
  try {
    *s = b->esvs.size();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_number_of_external_state_variables

mgis_status mgis_bv_get_external_state_variable_name(
    const char** n,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type i) {
  if (b == nullptr) {
    *n = nullptr;
    return mgis_report_failure("invalid argument");
  }
  try {
    const auto& ev = b->esvs.at(i);
    *n = ev.name.c_str();
    return mgis_report_success();
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
}  // end of mgis_bv_get_external_state_variable_name

void mgis_bv_free_behaviour(mgis_bv_Behaviour** d) {
  std::free(*d);
  *d = nullptr;
}  // end of mgis_bv_free_Behaviour

}  // end of extern "C"