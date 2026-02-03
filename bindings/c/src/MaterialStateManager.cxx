/*!
 * \file   bindings/c/src/MaterialStateManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/09/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/MaterialStateManager.h"

namespace mgis::behaviour::internals {

  [[nodiscard]] static MaterialStateManager::UpdatePolicy convert(
      const mgis_bv_MaterialStateManagerUpdatePolicy p) noexcept {
    if (p == MGIS_BV_UPDATE_POLICY) {
      return MaterialStateManager::UPDATE;
    }
    return MaterialStateManager::NOUPDATE;
  }  // end of convert

}  // end of namespace mgis::behaviour::internals

extern "C" {

mgis_status mgis_bv_material_state_manager_initializer_bind_gradients(
    mgis_bv_MaterialStateManagerInitializer* i,
    mgis_real* const p,
    mgis_size_type s) {
  if (i == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(material state manager initializer is null)");
  }
  if (p == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(memory pointer is null)");
  }
  try {
    i->gradients = std::span<mgis::real>(p, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_initializer_bind_gradients

mgis_status
mgis_bv_material_state_manager_initializer_bind_thermodynamic_forces(
    mgis_bv_MaterialStateManagerInitializer* i,
    mgis_real* const p,
    mgis_size_type s) {
  if (i == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(material state manager initializer is null)");
  }
  if (p == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(memory pointer is null)");
  }
  try {
    i->thermodynamic_forces = std::span<mgis::real>(p, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_initializer_bind_thermodynamic_forces

mgis_status
mgis_bv_material_state_manager_initializer_bind_internal_state_variables(
    mgis_bv_MaterialStateManagerInitializer* i,
    mgis_real* const p,
    mgis_size_type s) {
  if (i == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(material state manager initializer is null)");
  }
  if (p == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(memory pointer is null)");
  }
  try {
    i->internal_state_variables = std::span<mgis::real>(p, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_initializer_bind_internal_state_variables

mgis_status mgis_bv_material_state_manager_initializer_bind_stored_energies(
    mgis_bv_MaterialStateManagerInitializer* i,
    mgis_real* const p,
    mgis_size_type s) {
  if (i == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(material state manager initializer is null)");
  }
  if (p == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(memory pointer is null)");
  }
  try {
    i->stored_energies = std::span<mgis::real>(p, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_initializer_bind_stored_energies

mgis_status mgis_bv_material_state_manager_initializer_bind_dissipated_energies(
    mgis_bv_MaterialStateManagerInitializer* i,
    mgis_real* const p,
    mgis_size_type s) {
  if (i == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(material state manager initializer is null)");
  }
  if (p == nullptr) {
    return mgis_report_failure(
        "invalid argument "
        "(memory pointer is null)");
  }
  try {
    i->dissipated_energies = std::span<mgis::real>(p, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_initializer_bind_dissipated_energies

mgis_status mgis_bv_material_state_manager_get_number_of_integration_points(
    mgis_size_type* n, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *n = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *n = s->n;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_number_of_integration_points

mgis_status mgis_bv_material_state_manager_get_gradients(
    mgis_real** g, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *g = nullptr;
    return mgis_report_failure("null state manager");
  }
  *g = s->gradients.data();
  if (*g == nullptr) {
    return mgis_report_failure("no gradients defined");
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_gradients

mgis_status mgis_bv_material_state_manager_get_gradients_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->gradients_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_gradients_stride

mgis_status mgis_bv_material_state_manager_get_thermodynamic_forces(
    mgis_real** f, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *f = nullptr;
    return mgis_report_failure("null state manager");
  }
  *f = s->thermodynamic_forces.data();
  if (*f == nullptr) {
    return mgis_report_failure("no thermodynamic forces defined");
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_thermodynamic_forces

mgis_status mgis_bv_material_state_manager_get_thermodynamic_forces_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->thermodynamic_forces_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_thermodynamic_forces_stride

mgis_status mgis_bv_material_state_manager_get_internal_state_variables(
    mgis_real** ivs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *ivs = nullptr;
    return mgis_report_failure("null state manager");
  }
  *ivs = s->internal_state_variables.data();
  if (*ivs == nullptr) {
    return mgis_report_failure("no internal state variables defined");
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_internal_state_variables

mgis_status mgis_bv_material_state_manager_get_internal_state_variables_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->internal_state_variables_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_internal_state_variables_stride

mgis_status mgis_bv_material_state_manager_set_uniform_scalar_material_property(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    const mgis_real v) {
  return mgis_bv_material_state_manager_set_uniform_scalar_material_property2(
      m, n, v, MGIS_BV_UPDATE_POLICY);
}  // end of mgis_bv_material_state_manager_set_uniform_scalar_material_property

mgis_status
mgis_bv_material_state_manager_set_uniform_scalar_material_property2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    const mgis_real v,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    const auto b = mgis::behaviour::internals::convert(p);
    setMaterialProperty(*m, n, v, b);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_uniform_scalar_material_property

mgis_status mgis_bv_material_state_manager_set_uniform_material_property(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s) {
  return mgis_bv_material_state_manager_set_uniform_material_property2(
      m, n, v, s, MGIS_BV_UPDATE_POLICY);
}  // end of mgis_bv_material_state_manager_set_uniform_material_property

mgis_status mgis_bv_material_state_manager_set_uniform_material_property2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  if (v == nullptr) {
    return mgis_report_failure("null state values");
  }
  try {
    const auto& mp = getVariable(m->b.mps, n);
    const auto mpsize = getVariableSize(mp, m->b.hypothesis);
    const auto b = mgis::behaviour::internals::convert(p);
    if (s == MGIS_BV_LOCAL_STORAGE) {
      setMaterialProperty(*m, n, {v, static_cast<mgis_size_type>(mpsize)},
                          mgis::behaviour::MaterialStateManager::LOCAL_STORAGE,
                          b);
    } else {
      setMaterialProperty(
          *m, n, {v, static_cast<mgis_size_type>(mpsize)},
          mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE, b);
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_uniform_material_property2

mgis_status mgis_bv_material_state_manager_set_non_uniform_material_property(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s) {
  return mgis_bv_material_state_manager_set_non_uniform_material_property2(
      m, n, v, s, MGIS_BV_UPDATE_POLICY);
}  // end of mgis_bv_material_state_manager_set_non_uniform_material_property

mgis_status mgis_bv_material_state_manager_set_non_uniform_material_property2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  if (v == nullptr) {
    return mgis_report_failure("null material property values");
  }
  try {
    const auto& mp = getVariable(m->b.mps, n);
    const auto mpsize = getVariableSize(mp, m->b.hypothesis);
    const auto b = mgis::behaviour::internals::convert(p);
    if (s == MGIS_BV_LOCAL_STORAGE) {
      setMaterialProperty(
          *m, n, {v, static_cast<mgis_size_type>(m->n * mpsize)},
          mgis::behaviour::MaterialStateManager::LOCAL_STORAGE, b);
    } else {
      setMaterialProperty(
          *m, n, {v, static_cast<mgis_size_type>(m->n * mpsize)},
          mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE, b);
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_non_uniform_material_property

mgis_status mgis_bv_material_state_manager_is_material_property_defined(
    int* const b,
    const mgis_bv_MaterialStateManager* const m,
    const char* const n) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isMaterialPropertyDefined(*m, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_material_property_defined

mgis_status mgis_bv_material_state_manager_is_material_property_uniform(
    int* const b,
    const mgis_bv_MaterialStateManager* const m,
    const char* const n) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isMaterialPropertyUniform(*m, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_material_property_uniform

mgis_status mgis_bv_material_state_manager_set_uniform_scalar_mass_density(
    mgis_bv_MaterialStateManager* const m, const mgis_real v) {
  return mgis_bv_material_state_manager_set_uniform_scalar_mass_density2(
      m, v, MGIS_BV_UPDATE_POLICY);
}  // end of mgis_bv_material_state_manager_set_uniform_scalar_mass_density

mgis_status mgis_bv_material_state_manager_set_uniform_scalar_mass_density2(
    mgis_bv_MaterialStateManager* const m,
    const mgis_real v,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    const auto b = mgis::behaviour::internals::convert(p);
    setMassDensity(*m, v, b);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_uniform_scalar_mass_density

mgis_status mgis_bv_material_state_manager_set_non_uniform_mass_density(
    mgis_bv_MaterialStateManager* const m,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s) {
  return mgis_bv_material_state_manager_set_non_uniform_mass_density2(
      m, v, s, MGIS_BV_UPDATE_POLICY);
}

mgis_status mgis_bv_material_state_manager_set_non_uniform_mass_density2(
    mgis_bv_MaterialStateManager* const m,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  if (v == nullptr) {
    return mgis_report_failure("null state values");
  }
  try {
    const auto b = mgis::behaviour::internals::convert(p);
    if (s == MGIS_BV_LOCAL_STORAGE) {
      setMassDensity(*m, {v, static_cast<mgis_size_type>(m->n)},
                     mgis::behaviour::MaterialStateManager::LOCAL_STORAGE, b);
    } else {
      setMassDensity(*m, {v, static_cast<mgis_size_type>(m->n)},
                     mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE,
                     b);
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_non_uniform_mass_density

mgis_status mgis_bv_material_state_manager_is_mass_density_defined(
    int* const b, const mgis_bv_MaterialStateManager* const m) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isMassDensityDefined(*m);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_mass_density_defined

mgis_status mgis_bv_material_state_manager_is_mass_density_uniform(
    int* const b, const mgis_bv_MaterialStateManager* const m) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isMassDensityUniform(*m);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_mass_density_uniform

mgis_status
mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    const mgis_real v) {
  return mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable2(
      m, n, v, MGIS_BV_UPDATE_POLICY);
}  // end of
   // mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable

mgis_status
mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    const mgis_real v,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    const auto b = mgis::behaviour::internals::convert(p);
    setExternalStateVariable(*m, n, v, b);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable

mgis_status mgis_bv_material_state_manager_set_uniform_external_state_variable(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s) {
  return mgis_bv_material_state_manager_set_uniform_external_state_variable2(
      m, n, v, s, MGIS_BV_UPDATE_POLICY);
}  // end of mgis_bv_material_state_manager_set_uniform_external_state_variable

mgis_status mgis_bv_material_state_manager_set_uniform_external_state_variable2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  if (v == nullptr) {
    return mgis_report_failure("null state values");
  }
  try {
    const auto& esv = getVariable(m->b.esvs, n);
    const auto esvsize = getVariableSize(esv, m->b.hypothesis);
    const auto b = mgis::behaviour::internals::convert(p);
    if (s == MGIS_BV_LOCAL_STORAGE) {
      setExternalStateVariable(
          *m, n, {v, static_cast<mgis_size_type>(esvsize)},
          mgis::behaviour::MaterialStateManager::LOCAL_STORAGE, b);
    } else {
      setExternalStateVariable(
          *m, n, {v, static_cast<mgis_size_type>(esvsize)},
          mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE, b);
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_set_uniform_external_state_variable

mgis_status
mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s) {
  return mgis_bv_material_state_manager_set_non_uniform_external_state_variable2(
      m, n, v, s, MGIS_BV_UPDATE_POLICY);
}

mgis_status
mgis_bv_material_state_manager_set_non_uniform_external_state_variable2(
    mgis_bv_MaterialStateManager* const m,
    const char* const n,
    mgis_real* const v,
    const mgis_bv_MaterialStateManagerStorageMode s,
    const mgis_bv_MaterialStateManagerUpdatePolicy p) {
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  if (v == nullptr) {
    return mgis_report_failure("null state values");
  }
  try {
    const auto& esv = getVariable(m->b.esvs, n);
    const auto esvsize = getVariableSize(esv, m->b.hypothesis);
    const auto b = mgis::behaviour::internals::convert(p);
    if (s == MGIS_BV_LOCAL_STORAGE) {
      setExternalStateVariable(
          *m, n, {v, static_cast<mgis_size_type>(esvsize * m->n)},
          mgis::behaviour::MaterialStateManager::LOCAL_STORAGE, b);
    } else {
      setExternalStateVariable(
          *m, n, {v, static_cast<mgis_size_type>(esvsize * m->n)},
          mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE, b);
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of
   // mgis_bv_material_state_manager_set_non_uniform_external_state_variable

mgis_status mgis_bv_material_state_manager_is_external_state_variable_defined(
    int* const b,
    const mgis_bv_MaterialStateManager* const m,
    const char* const n) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isExternalStateVariableDefined(*m, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_external_state_variable_defined

mgis_status mgis_bv_material_state_manager_is_external_state_variable_uniform(
    int* const b,
    const mgis_bv_MaterialStateManager* const m,
    const char* const n) {
  *b = 0;
  if (m == nullptr) {
    return mgis_report_failure("null state manager");
  }
  try {
    *b = isExternalStateVariableUniform(*m, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_is_external_state_variable_uniform

}  // end of extern "C"
