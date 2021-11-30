/*!
 * \file   bindings/c/include/MGIS/Behaviour/MaterialStateManager.h
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

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_H
#define LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_H

#include "MGIS/Config.h"
#include "MGIS/Behaviour/Behaviour.h"

#ifdef __cplusplus
#include "MGIS/Behaviour/MaterialStateManager.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_MaterialStateManagerInitializer =
    mgis::behaviour::MaterialStateManagerInitializer;
using mgis_bv_MaterialStateManager = mgis::behaviour::MaterialStateManager;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_MaterialStateManagerInitializer
    mgis_bv_MaterialStateManagerInitializer;
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_MaterialStateManager mgis_bv_MaterialStateManager;
#endif

typedef enum {
  MGIS_BV_LOCAL_STORAGE = 0,
  MGIS_BV_EXTERNAL_STORAGE = 1
} mgis_bv_MaterialStateManagerStorageMode;
/*!
 * \brief bind the gradients to the given array
 * \param[in,out] s: initializer
 * \param[in] g: pointer to a memory area meant to store the gradients values
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status mgis_bv_material_state_manager_initializer_bind_gradients(
    mgis_bv_MaterialStateManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \brief bind the thermodynamic forces to the given array
 * \param[in,out] s: initializer
 * \param[in] p: pointer to a memory area meant to store the thermodynamic
 * forces values
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status
mgis_bv_material_state_manager_initializer_bind_thermodynamic_forces(
    mgis_bv_MaterialStateManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \brief bind the internal state variables to the given array
 * \param[in,out] s: initializer
 * \param[in] p: pointer to a memory area meant to store the internal
 * state variables values
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status
mgis_bv_material_state_manager_initializer_bind_internal_state_variables(
    mgis_bv_MaterialStateManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \brief bind the stored energies to the given array
 * \param[in,out] s: initializer
 * \param[in] p: pointer to a memory area meant to store the stored energies
 * values.
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status mgis_bv_material_state_manager_initializer_bind_stored_energies(
    mgis_bv_MaterialStateManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \brief bind the dissipated energies to the given array
 * \param[in,out] s: initializer
 * \param[in] p: pointer to a memory area meant to store the dissipated energies
 * values.
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status mgis_bv_material_state_manager_initializer_bind_dissipated_energies(
    mgis_bv_MaterialStateManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \param[out] n: number of integration points
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_number_of_integration_points(
    mgis_size_type*, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] g: a pointer to the array of gradients
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_material_state_manager_get_gradients(
    mgis_real**, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] gs: a pointer where the gradients' stride will be written
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_material_state_manager_get_gradients_stride(
    mgis_size_type* const, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] f: a pointer to the array of thermodynamic_forces
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_thermodynamic_forces(
    mgis_real**, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] fs: a pointer where the thermodynamic forces' stride will be
 * written.
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_thermodynamic_forces_stride(
    mgis_size_type* const, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] ivs: a pointer to the array of internal state variables
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_internal_state_variables(
    mgis_real**, mgis_bv_MaterialStateManager* const);
/*!
 * \param[out] fs: a pointer where the internal state variables' stride will be
 * written.
 * \param[in]  s: state manager
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_internal_state_variables_stride(
    mgis_size_type* const, mgis_bv_MaterialStateManager* const);
/*!
 * \brief set the value of an uniform scalar material property
 * \param[in] s: state manager
 * \param[in] n: name of the material property
 * \param[in] v: value of the material property
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_uniform_scalar_material_property(
    mgis_bv_MaterialStateManager* const, const char* const, const mgis_real);
/*!
 * \brief set the value of an uniform material property
 * \param[in] s: state manager
 * \param[in] n: name of the material property
 * \param[in] v: values of the material property
 * \param[in] sm: storage mode
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_uniform_material_property(
    mgis_bv_MaterialStateManager* const,
    const char* const,
    mgis_real* const,
    const mgis_bv_MaterialStateManagerStorageMode);
/*!
 * \brief set the value of an uniform material property
 * \param[in] s: state manager
 * \param[in] n: name of the material property
 * \param[in] v: values of the material property
 * \param[in] sm: storage mode
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_non_uniform_material_property(
    mgis_bv_MaterialStateManager* const,
    const char* const,
    mgis_real* const,
    const mgis_bv_MaterialStateManagerStorageMode);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_is_material_property_defined(
    int* const, const mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_is_material_property_uniform(
    int* const, const mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_uniform_material_property(
    mgis_real* const, mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_non_uniform_material_property(
    mgis_real** const, mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
    mgis_bv_MaterialStateManager* const, const char* const, const mgis_real);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_uniform_external_state_variable(
    mgis_bv_MaterialStateManager* const,
    const char* const,
    mgis_real* const,
    const mgis_bv_MaterialStateManagerStorageMode);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
    mgis_bv_MaterialStateManager* const,
    const char* const,
    mgis_real* const,
    const mgis_bv_MaterialStateManagerStorageMode);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_is_external_state_variable_defined(
    int* const, const mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_is_external_state_variable_uniform(
    int* const, const mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_uniform_external_state_variable(
    mgis_real* const, mgis_bv_MaterialStateManager* const, const char* const);

MGIS_C_EXPORT mgis_status
mgis_bv_material_state_manager_get_non_uniform_external_state_variable(
    mgis_real** const, mgis_bv_MaterialStateManager* const, const char* const);

#ifdef __cplusplus
}  // end of extern "C"
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_H */
