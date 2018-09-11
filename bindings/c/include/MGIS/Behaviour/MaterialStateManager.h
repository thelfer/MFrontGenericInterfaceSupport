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
using mgis_bv_MaterialStateManager = mgis::behaviour::MaterialStateManager;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_MaterialStateManager mgis_bv_MaterialStateManager;
#endif

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
 * \param[out] f: a pointer to the array of internal state variables
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

#ifdef __cplusplus
}  // end of extern "C"
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_H */
