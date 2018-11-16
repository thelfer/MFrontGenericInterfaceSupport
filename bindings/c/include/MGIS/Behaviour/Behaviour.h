/*!
 * \file   bindings/c/include/MFront/Behaviour/Behaviour.h
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

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOUR_H
#define LIB_MGIS_BEHAVIOUR_BEHAVIOUR_H

#ifdef __cplusplus
#include <cstddef>
#else /*  __cplusplus */
#include <stddef.h>
#endif /*  __cplusplus */

#include "MGIS/Config.h"
#include "MGIS/Status.h"
#include "MGIS/Behaviour/Variable.h"

#ifdef __cplusplus
#include "MGIS/Behaviour/Behaviour.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_Behaviour = mgis::behaviour::Behaviour;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_Behaviour mgis_bv_Behaviour;
#endif

typedef enum {
  MGIS_BV_ISOTROPIC = 0,
  MGIS_BV_ORTHOTROPIC = 1
} mgis_bv_BehaviourSymmetry;

typedef enum {
  MGIS_BV_GENERALBEHAVIOUR = 0,
  MGIS_BV_STANDARDSTRAINBASEDBEHAVIOUR = 1,
  MGIS_BV_STANDARDFINITESTRAINBEHAVIOUR = 2,
  MGIS_BV_COHESIVEZONEMODEL = 3
} mgis_bv_BehaviourType;

//! kinematic of the behaviour treated
typedef enum {
  MGIS_BV_UNDEFINEDKINEMATIC = 0,
  MGIS_BV_SMALLSTRAINKINEMATIC = 1,
  MGIS_BV_COHESIVEZONEKINEMATIC = 2,
  MGIS_BV_FINITESTRAINKINEMATIC_F_CAUCHY = 3,
  MGIS_BV_FINITESTRAINKINEMATIC_ETO_PK1 = 4
} mgis_bv_BehaviourKinematic;

/*!
 * \brief load a behaviour
 * \param[out] ptr: behaviour
 * \param[in] l: library name
 * \param[in] b: behaviour name
 * \param[in] h: hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_bv_load_behaviour(mgis_bv_Behaviour**,
                                                 const char* const,
                                                 const char* const,
                                                 const char* const);
/*!
 * \brief retrieve the library
 * \param[out] l: library
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_library(const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the source
 * \param[out] s: source
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_source(const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the hypothesis
 * \param[out] h: hypothesis
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_hypothesis(const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the behaviour name
 * \param[out] b: behaviour name
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_name(
    const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the function name
 * \param[out] b: function name
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_function_name(
    const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the `TFEL` version used to generate the behaviour
 * \param[out] v: version
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_tfel_version(
    const char**, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the behaviour symmetry
 * \param[out] s: symmetry
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_symmetry(
    mgis_bv_BehaviourSymmetry* const, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the behaviour type
 * \param[out] t: behaviour type
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_type(
    mgis_bv_BehaviourType* const, const mgis_bv_Behaviour* const);
/*!
 * \brief retrieve the behaviour kinematic
 * \param[out] k: behaviour kinematic
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_behaviour_kinematic(
    mgis_bv_BehaviourKinematic* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the number of material properties
 * \param[out] c: number of the material properties
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_number_of_material_properties(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of a material property
 * \param[out] c: material property name
 * \param[in] b: behaviour
 * \param[in] i: material property index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_material_property_name(
    const char**, const mgis_bv_Behaviour* const, const mgis_size_type);
/*!
 * \brief return the type of an material property
 * \param[out] t: material property type
 * \param[in] b: behaviour
 * \param[in] i: material property index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_material_property_type(
    mgis_bv_VariableType* const,
    const mgis_bv_Behaviour* const,
    const mgis_size_type);
/*!
 * \brief return the number of internal state variables
 * \param[out] c: number of the internal state variables
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_number_of_internal_state_variables(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of an internal state variable
 * \param[out] c: internal state variable name
 * \param[in] b: behaviour
 * \param[in] i: internal state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_internal_state_variable_name(
    const char**, const mgis_bv_Behaviour* const, const mgis_size_type);
/*!
 * \brief return the name of an internal state variable
 * \param[out] o: internal state variable offset
 * \param[in]  b: behaviour
 * \param[in]  n: internal state variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_internal_state_variable_offset(
    mgis_size_type* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return the type of an internal state variable
 * \param[out] t: internal state variable type
 * \param[in] b: behaviour
 * \param[in] i: internal state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_internal_state_variable_type(
    mgis_bv_VariableType* const,
    const mgis_bv_Behaviour* const,
    const mgis_size_type);
/*!
 * \brief return the number of external state variables
 * \param[out] c: number of the external state variables
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_number_of_external_state_variables(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of an external state variable
 * \param[out] c: external state variable name
 * \param[in] b: behaviour
 * \param[in] i: external state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_external_state_variable_name(
    const char **, const mgis_bv_Behaviour *const, const mgis_size_type);
/*!
 * \brief return the type of an external state variable
 * \param[out] t: external state variable type
 * \param[in] b: behaviour
 * \param[in] i: external state variable index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_external_state_variable_type(
    mgis_bv_VariableType* const,
    const mgis_bv_Behaviour* const,
    const mgis_size_type);
/*!
 * \brief return the number of parameters
 * \param[out] c: number of the parameters
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_number_of_parameters(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of an parameter
 * \param[out] c: parameter name
 * \param[in] b: behaviour
 * \param[in] i: parameter index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_parameter_name(
    const char**, const mgis_bv_Behaviour* const, const mgis_size_type);
/*!
 * \brief set the value of a parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 * \param[in] v: parameter value
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_set_parameter(
    const mgis_bv_Behaviour *const, const char *const, const double);
/*!
 * \brief return the number of integer parameters
 * \param[out] c: number of the integer parameters
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_number_of_integer_parameters(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of an integer parameter
 * \param[out] c: integer parameter name
 * \param[in] b: behaviour
 * \param[in] i: integer parameter index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_integer_parameter_name(
    const char**, const mgis_bv_Behaviour* const, const mgis_size_type);
/*!
 * \brief set the value of a parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 * \param[in] v: parameter value
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_set_integer_parameter(
    const mgis_bv_Behaviour *const, const char *const, const int);
/*!
 * \brief return the number of unsigned short parameters
 * \param[out] c: number of the unsigned short parameters
 * \param[in] b: behaviour
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_number_of_unsigned_short_parameters(
    mgis_size_type* const, const mgis_bv_Behaviour* const);
/*!
 * \brief return the name of an unsigned short parameter
 * \param[out] c: unsigned short parameter name
 * \param[in] b: behaviour
 * \param[in] i: unsigned short parameter index
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_unsigned_short_parameter_name(
    const char**, const mgis_bv_Behaviour* const, const mgis_size_type);
/*!
 * \brief set the value of a parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 * \param[in] v: parameter value
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_set_unsigned_short_parameter(
    const mgis_bv_Behaviour *const, const char *const, const unsigned short);
/*!
 * \brief grant access to the default value of a parameter
 * \param[out] v: default value of the parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_parameter_default_value(
    double *const, const mgis_bv_Behaviour *const, const char *const);
/*!
 * \brief grant access to the default value of a parameter
 * \param[out] v: default value of the parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_integer_parameter_default_value(
    int *const, const mgis_bv_Behaviour *const, const char *const);
/*!
 * \brief grant access to the default value of a parameter
 * \param[out] v: default value of the parameter
 * \param[in] b: behaviour description
 * \param[in] n: parameter name
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_get_unsigned_short_parameter_default_value(
    unsigned short *const, const mgis_bv_Behaviour *const, const char *const);
/*!
 * \brief return true if the given variable has bounds
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_bounds(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return true if the given variable has a lower bound
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_lower_bound(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return true if the given variable has a upper bound
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_upper_bound(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return the lower bound of the given variable
 * \param[out] v: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_lower_bound(
    long double* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return the upper bound of the given variable
 * \param[out] v: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_upper_bound(
    long double* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return true if the given variable has physical bounds
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_physical_bounds(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return true if the given variable has a lower physical bound
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_lower_physical_bound(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return true if the given variable has a upper physical bound
 * \param[out] r: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_has_upper_physical_bound(
    int* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return the lower physical bound of the given variable
 * \param[out] v: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_lower_physical_bound(
    long double* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief return the upper physical bound of the given variable
 * \param[out] v: returned value
 * \param[in] b: behaviour
 * \param[in] n: variable name
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_get_upper_physical_bound(
    long double* const, const mgis_bv_Behaviour* const, const char* const);
/*!
 * \brief free the memory associated with the given behaviour.
 * \param[in,out] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_free_behaviour(mgis_bv_Behaviour**);

#ifdef __cplusplus
}  // end of extern "C"
#endif

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOUR_H */
