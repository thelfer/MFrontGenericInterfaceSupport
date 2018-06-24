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

MGIS_C_EXPORT mgis_status
mgis_behaviour_load_description(mgis_behaviour_Description**,
                                const char* const,
                                const char* const,
                                const char* const);

MGIS_C_EXPORT void mgis_behaviour_free_description(
    mgis_behaviour_Description* const);

/*!
 * \brief return the numer of material properties
 * \param[in] b: behaviour description
 */
MGIS_C_EXPORT mgis_status mgis_behaviour_getNumerOfMaterialProperties(
    mgis_size_type* const, const mgis_behaviour_Description* const);

#ifdef __cplusplus
}
#endif

#endif /* LIB_MGIS_BEHAVIOUR_DESCRIPTION_H */
