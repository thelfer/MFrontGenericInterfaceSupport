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

#ifdef __cplusplus
#include "MGIS/Behaviour/State.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_State = mgis::behaviour::State;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_State mgis_bv_State;
#endif

#ifdef __cplusplus
}  // end of extern "C"
#endif

#endif /* LIB_MGIS_BEHAVIOUR_STATE_H */
