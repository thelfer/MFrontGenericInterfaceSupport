/*!
 * \file   bindings/c/include/MFront/Behaviour/Description.h
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

#ifndef LIB_MGIS_BEHAVIOUR_DESCRIPTION_H
#define LIB_MGIS_BEHAVIOUR_DESCRIPTION_H

#ifdef __cplusplus
#include <cstddef>
#include "MGIS/Behaviour/Description.hxx"
extern "C" {
#else /*  __cplusplus */
#include <stddef.h>
#endif /*  __cplusplus */

#include "MGIS/Config-c.h"
#include "MGIS/Status.h"
#include "MGIS/Behaviour/Variable.h"

#ifdef __cplusplus
//! a simple alias
using mgis_bv_Description = mgis::behaviour::Description;
#else
/*!
 * \brief an opaque structure which can only be accessed through the mgis API.
 */
typedef struct mgis_bv_Description mgis_bv_Description;
#endif

typedef enum {
  MGIS_BV_ISOTROPIC,
  MGIS_BV_ORTHOTROPIC
} mgis_bv_BehaviourSymmetry;

typedef enum {
  MGIS_BV_GENERALBEHAVIOUR,
  MGIS_BV_STANDARDSTRAINBASEDBEHAVIOUR,
  MGIS_BV_STANDARDFINITESTRAINBEHAVIOUR,
  MGIS_BV_COHESIVEZONEMODEL
} mgis_bv_BehaviourType;

//! kinematic of the behaviour treated
typedef enum {
  MGIS_BV_UNDEFINEDKINEMATIC,
  MGIS_BV_SMALLSTRAINKINEMATIC,
  MGIS_BV_COHESIVEZONEKINEMATIC,
  MGIS_BV_FINITESTRAINKINEMATIC_F_CAUCHY,
  MGIS_BV_FINITESTRAINKINEMATIC_ETO_PK1
} mgis_bv_BehaviourKinematic;

/*!
 * \brief load a behaviour description
 * \param[out] d: behaviour description
 * \param[in] l: library name
 * \param[in] b: behaviour name
 * \param[in] h: hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_bv_load_description(mgis_bv_Description**,
                                                   const char* const,
                                                   const char* const,
                                                   const char* const);
/*!
 * \brief retrieve the library
 * \param[out] l: library
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_library(const char**,
                                              const mgis_bv_Description* const);
/*!
 * \brief retrieve the source
 * \param[out] s: source
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_source(const char**,
                                             const mgis_bv_Description* const);
/*!
 * \brief retrieve the hypothesis
 * \param[out] h: hypothesis
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status
mgis_bv_get_hypothesis(const char**, const mgis_bv_Description* const);
/*!
 * \brief retrieve the behaviour name
 * \param[out] b: behaviour name
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status
mgis_bv_get_behaviour_name(const char**, const mgis_bv_Description* const);
/*!
 * \brief retrieve the function name
 * \param[out] b: function name
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status
mgis_bv_get_function_name(const char**, const mgis_bv_Description* const);
/*!
 * \brief retrieve the `TFEL` version used to generate the behaviour
 * \param[out] v: version
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status
mgis_bv_get_tfel_version(const char**, const mgis_bv_Description* const);
/*!
 * \brief retrieve the behaviour symmetry
 * \param[out] s: symmetry
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_symmetry(
    mgis_bv_BehaviourSymmetry* const, const mgis_bv_Description* const);
/*!
 * \brief retrieve the behaviour type
 * \param[out] t: behaviour type
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_type(
    mgis_bv_BehaviourType* const, const mgis_bv_Description* const);
/*!
 * \brief retrieve the behaviour kinematic
 * \param[out] k: behaviour kinematic
 * \param[in] d: description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_behaviour_kinematic(
    mgis_bv_BehaviourKinematic* const, const mgis_bv_Description* const);
/*!
 * \brief return the number of material properties
 * \param[out] c: number of the material properties
 * \param[in] b: behaviour description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_number_of_material_properties(
    mgis_size_type* const, const mgis_bv_Description* const);
/*!
 * \brief return the name of a material property
 * \param[out] c: material property name
 * \param[in] b: behaviour description
 * \param[in] i: material property index
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_material_property_name(
    const char**, const mgis_bv_Description* const, const mgis_size_type);
/*!
 * \brief return the number of internal state variables
 * \param[out] c: number of the internal state variables
 * \param[in] b: behaviour description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_number_of_internal_state_variables(
    mgis_size_type* const, const mgis_bv_Description* const);
/*!
 * \brief return the name of an internal state variable
 * \param[out] c: internal state variable name
 * \param[in] b: behaviour description
 * \param[in] i: internal state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_internal_state_variable_name(
    const char**, const mgis_bv_Description* const, const mgis_size_type);
/*!
 * \brief return the type of an internal state variable
 * \param[out] t: internal state variable type
 * \param[in] b: behaviour description
 * \param[in] i: internal state variable index
 */
MGIS_C_EXPORT mgis_status
mgis_bv_get_internal_state_variable_type(mgis_bv_VariableType* const,
                                         const mgis_bv_Description* const,
                                         const mgis_size_type);

/*!
 * \brief return the number of external state variables
 * \param[out] c: number of the external state variables
 * \param[in] b: behaviour description
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_number_of_external_state_variables(
    mgis_size_type* const, const mgis_bv_Description* const);
/*!
 * \brief return the name of an external state variable
 * \param[out] c: external state variable name
 * \param[in] b: behaviour description
 * \param[in] i: external state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_external_state_variable_name(
    const char**, const mgis_bv_Description* const, const mgis_size_type);
/*!
 * \brief free the memory associated with the given description.
 * \param[in,out] d: description
 */
MGIS_C_EXPORT void mgis_bv_free_description(mgis_bv_Description**);

#ifdef __cplusplus
}
#endif

#endif /* LIB_MGIS_BEHAVIOUR_DESCRIPTION_H */
