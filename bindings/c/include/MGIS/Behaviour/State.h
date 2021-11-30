/*!
 * \file   include/MGIS/Behaviour/State.h
 * \brief
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_STATE_H
#define LIB_MGIS_BEHAVIOUR_STATE_H

#include "MGIS/Status.h"

#ifdef __cplusplus
#include "MGIS/Behaviour/State.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_State = mgis::behaviour::State;
#else /*  __cplusplus */
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_State mgis_bv_State;
#endif /*  __cplusplus */

/*!
 * \brief set a gradient' value in a state
 * \param[out] s: state
 * \param[in]  n: name
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_gradient_by_name(mgis_bv_State* const,
                                   const char* const,
                                   const mgis_real* const);
/*!
 * \brief get a gradient' value in a state
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_get_gradient_by_name(mgis_real**,
                                   mgis_bv_State* const,
                                   const char* const);
/*!
 * \brief set a gradient' value in a state
 * \param[out] s: state
 * \param[in]  o: offset
 * \param[in]  n: variable size
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_gradient_by_offset(mgis_bv_State* const,
                                     const mgis_size_type,
                                     const mgis_size_type,
                                     const mgis_real* const);
/*!
 * \brief get a gradient' value in a state
 * \param[out] v: pointer to gradient' value(s)
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_gradient_by_offset(
    mgis_real**, mgis_bv_State* const, const mgis_size_type);
/*!
 * \brief set a thermodynamic force' value in a state
 * \param[out] s: state
 * \param[in]  n: name
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_thermodynamic_force_by_name(mgis_bv_State* const,
                                              const char* const,
                                              const mgis_real* const);
/*!
 * \brief get a thermodynamic force' value in a state
 * \param[out] v: pointer to thermodynamic force' value(s)
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_get_thermodynamic_force_by_name(mgis_real**,
                                              mgis_bv_State* const,
                                                 const char* const);
/*!
 * \brief set a thermodynamic force' value in a state
 * \param[out] s: state
 * \param[in]  o: offset
 * \param[in]  n: variable size
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_thermodynamic_force_by_offset(mgis_bv_State* const,
                                                const mgis_size_type,
                                                const mgis_size_type,
                                                const mgis_real* const);
/*!
 * \brief get a thermodynamic force' value in a state
 * \param[out] v: pointer to thermodynamic force' value(s)
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_thermodynamic_force_by_offset(
    mgis_real**, mgis_bv_State* const, const mgis_size_type);
/*!
 * \brief set a material property' value in a state
 * \param[out] s: state
 * \param[in]  n: name
 * \param[in]  v: value
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_material_property_by_name(mgis_bv_State* const,
                                            const char* const,
                                            const mgis_real);
/*!
 * \brief get a material property' value in a state
 * \param[out] v: material property' value
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_get_material_property_by_name(mgis_real** const,
                                            mgis_bv_State* const,
                                            const char* const);
/*!
 * \brief set a material property' value in a state
 * \param[out] s: state
 * \param[in] o: offset
 * \param[in] v: value
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_set_material_property_by_offset(
    mgis_bv_State* const, const mgis_size_type, const mgis_real);
/*!
 * \brief get a material property' value in a state
 * \param[out] v: material property' value
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_material_property_by_offset(
    mgis_real** const, mgis_bv_State* const, const mgis_size_type);
/*!
 * \brief set a internal state variable' value in a state
 * \param[out] s: state
 * \param[in]  n: name
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_set_internal_state_variable_by_name(
    mgis_bv_State* const,
    const char* const,
    const mgis_real* const);
/*!
 * \brief get a pointer to the state variables
 * \param[out] v: pointer to internal state variables
 * \param[in] s: state
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_internal_state_variables(
    mgis_real**,
    mgis_bv_State* const);
/*!
 * \brief get a internal state variable' value in a state
 * \param[out] v: pointer to internal state variable' value(s)
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_internal_state_variable_by_name(
    mgis_real**,
    mgis_bv_State* const,
    const char* const);
/*!
 * \brief set a internal state variable' value in a state
 * \param[out] s: state
 * \param[in]  o: offset
 * \param[in]  n: variable size
 * \param[in]  v: value(s)
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_internal_state_variable_by_offset(mgis_bv_State* const,
                                                    const mgis_size_type,
                                                    const mgis_size_type,
                                                    const mgis_real* const);
/*!
 * \brief get a internal state variable' value in a state
 * \param[out] v: pointer to internal state variable' value(s)
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_internal_state_variable_by_offset(
    mgis_real**, mgis_bv_State* const, const mgis_size_type);
/*!
 * \brief set a external state variable' value in a state
 * \param[out] s: state
 * \param[in] n: name
 * \param[in] v: value
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_scalar_external_state_variable_by_name(mgis_bv_State* const,
                                                         const char* const,
                                                         const mgis_real);
/*!
 * \brief set a external state variable' value in a state
 * \param[out] s: state
 * \param[in] n: name
 * \param[in] v: value
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_set_external_state_variable_by_name(
    mgis_bv_State* const, const char* const, const mgis_real* const);
/*!
 * \brief get a pointer to the state variables
 * \param[out] v: pointer to external state variables
 * \param[in] s: state
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_external_state_variables(
    mgis_real**,
    mgis_bv_State* const);
/*!
 * \brief get a external state variable' value in a state
 * \param[out] v: external state variable' value
 * \param[in] s: state
 * \param[in] n: name
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_external_state_variable_by_name(
    mgis_real** const,
    mgis_bv_State* const,
    const char* const);
/*!
 * \brief set a external state variable' value in a state
 * \param[out] s: state
 * \param[in] o: offset
 * \param[in] v: value
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_scalar_external_state_variable_by_offset(mgis_bv_State* const,
                                                           const mgis_size_type,
                                                           const mgis_real);
/*!
 * \brief set a external state variable' value in a state
 * \param[out] s: state
 * \param[in] o: offset
 * \param[in] v: value
 * \param[in] vs: number of values
 */
MGIS_C_EXPORT mgis_status
mgis_bv_state_set_external_state_variable_by_offset(mgis_bv_State* const,
                                                    const mgis_size_type,
                                                    const mgis_real* const,
                                                    const mgis_size_type);
/*!
 * \brief get a external state variable' value in a state
 * \param[out] v: external state variable' value
 * \param[in]  s: state
 * \param[in]  o: offset
 */
MGIS_C_EXPORT mgis_status mgis_bv_state_get_external_state_variable_by_offset(
    mgis_real** const, mgis_bv_State* const, const mgis_size_type);

#ifdef __cplusplus
}  // end of extern "C"
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_STATE_H */
