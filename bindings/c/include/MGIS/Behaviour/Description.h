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
#else  /*  __cplusplus */
#include <stddef.h>
#endif /*  __cplusplus */

#include "MGIS/Config-c.h"
#include "MGIS/Status.h"

#ifdef __cplusplus
//! a simple alias
using mgis_behaviour_Description = mgis::behaviour::Description;
#else
/*!
 * \brief an opaque structure which can only be accessed through the mgis API.
 */
typedef struct mgis_behaviour_Description mgis_behaviour_Description;
#endif

/*!
 * \brief load a behaviour description
 * \param[out] d: behaviour description
 * \param[in] l: library name
 * \param[in] b: behaviour name
 * \param[in] h: hypothesis
 */
MGIS_C_EXPORT mgis_status
mgis_behaviour_load_description(mgis_behaviour_Description**,
                                const char* const,
                                const char* const,
                                const char* const);
/*!
 * \brief return the number of material properties
 * \param[out] c: number of the material properties
 * \param[in] b: behaviour description
 */
MGIS_C_EXPORT mgis_status mgis_behaviour_get_number_of_material_properties(
    mgis_size_type *const, const mgis_behaviour_Description *const);
/*!
 * \brief return the numer of material properties
 * \param[in] c: material property name
 * \param[in] b: behaviour description
 * \param[in] i: material property index
 */
MGIS_C_EXPORT mgis_status mgis_behaviour_get_material_property_name(
    const char **, const mgis_behaviour_Description *const,
    const mgis_size_type);
/*!
 * \brief free the memory associated with the given description.
 * \param[in,out] d: description
 */
MGIS_C_EXPORT void
mgis_behaviour_free_description(mgis_behaviour_Description **);

#ifdef __cplusplus
}
#endif

#endif /* LIB_MGIS_BEHAVIOUR_DESCRIPTION_H */
