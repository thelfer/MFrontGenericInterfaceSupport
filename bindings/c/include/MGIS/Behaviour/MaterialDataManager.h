/*!
 * \file   MaterialDataManager.h
 * \brief
 * \author Thomas Helfer
 * \date   05/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_H
#define LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_H

#include "MGIS/Config.h"
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/MaterialStateManager.h"

#ifdef __cplusplus
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_MaterialDataManagerInitializer =
    mgis::behaviour::MaterialDataManagerInitializer;
using mgis_bv_MaterialDataManager = mgis::behaviour::MaterialDataManager;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_MaterialDataManagerInitializer
    mgis_bv_MaterialDataManagerInitializer;
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_MaterialDataManager mgis_bv_MaterialDataManager;
#endif

/*!
 * \param[out] d: a pointer to the created data manager initializer
 */
MGIS_C_EXPORT mgis_status mgis_bv_create_material_data_manager_initializer(
    mgis_bv_MaterialDataManagerInitializer**);
/*!
 * \brief bind the tangent operator to the given array
 * \param[in,out] d: initializer
 * \param[in] K: pointer to a memory area meant to store the consistent tangent operator values
 * \param[in] s: size of the memory area
 */
MGIS_C_EXPORT
mgis_status mgis_bv_material_data_manager_initializer_bind_tangent_operator(
    mgis_bv_MaterialDataManagerInitializer*, mgis_real* const, mgis_size_type);
/*!
 * \brief set the state at the beginning of the time step
 * \param[out] s: pointer to a pointer to the state initializer
 * \param[in]  d: material data manager initalizer
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_data_manager_initializer_get_state_0_initializer(
    mgis_bv_MaterialStateManagerInitializer**,
    mgis_bv_MaterialDataManagerInitializer* const);
/*!
 * \brief set the state at the beginning of the time step
 * \param[out] s: pointer to a pointer to the state initializer
 * \param[in]  d: material data manager initalizer
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_data_manager_initializer_get_state_1_initializer(
    mgis_bv_MaterialStateManagerInitializer**,
    mgis_bv_MaterialDataManagerInitializer* const);
/*!
 * \brief free the memory associated with the given material data manager
 * initializer.
 * \param[in,out] d: data manager initializer
 */
MGIS_C_EXPORT mgis_status mgis_bv_free_material_data_manager_initializer(
    mgis_bv_MaterialDataManagerInitializer**);

/*!
 * \param[out] d: a pointer to the created data manager
 * \param[in]  b: behaviour
 * \param[in]  n: number of integration points
 */
MGIS_C_EXPORT mgis_status mgis_bv_create_material_data_manager(
    mgis_bv_MaterialDataManager**, const mgis_bv_Behaviour* const,
    const mgis_size_type);
/*!
 * \param[out] d: a pointer to the created data manager
 * \param[in]  b: behaviour
 * \param[in]  n: number of integration points
 * \param[in]  i: initializer
 */
MGIS_C_EXPORT mgis_status mgis_bv_create_material_data_manager_with_initializer(
    mgis_bv_MaterialDataManager**,
    const mgis_bv_Behaviour* const,
    const mgis_size_type,
    const mgis_bv_MaterialDataManagerInitializer* const);
/*!
 * \brief set the state at the beginning of the time step
 * \param[out] s: pointer to a pointer to the state
 * \param[in]  d: material data manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_material_data_manager_get_state_0(
    mgis_bv_MaterialStateManager**, mgis_bv_MaterialDataManager * const);
/*!
 * \brief set the state at the end of the time step
 * \param[out] s: pointer to a pointer to the state
 * \param[in]  d: material data manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_material_data_manager_get_state_1(
    mgis_bv_MaterialStateManager**, mgis_bv_MaterialDataManager * const);
/*!
 * \brief return the time step scaling factor
 * \param[out] rdt: time step scaling factor
 * \param[in]  d: behaviour data
 */
MGIS_C_EXPORT mgis_status
mgis_bv_material_data_manager_get_time_step_scaling_factor(
    mgis_real * const, const mgis_bv_MaterialDataManager* const);
/*!
 * \brief return the tangent operator
 * \param[out] K: tangent operator
 * \param[in]  d: behaviour data
 */
MGIS_C_EXPORT mgis_status mgis_bv_material_data_manager_get_tangent_operator(
    mgis_real * *const, mgis_bv_MaterialDataManager* const);
/*!
 * \brief update the material data manager
 * \param[in,out] d: data manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_update_material_data_manager(
    mgis_bv_MaterialDataManager * const);
/*!
 * \brief revert the material data manager
 * \param[in,out] d: data manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_revert_material_data_manager(
    mgis_bv_MaterialDataManager * const);
/*!
 * \brief free the memory associated with the given material data manager.
 * \param[in,out] d: data manager
 */
MGIS_C_EXPORT mgis_status mgis_bv_free_material_data_manager(
    mgis_bv_MaterialDataManager**);

#ifdef __cplusplus
}  // end of extern "C"
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_H */
