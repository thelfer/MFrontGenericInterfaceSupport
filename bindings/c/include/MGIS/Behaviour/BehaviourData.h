/*!
 * \file   include/MGIS/Behaviour/BehaviourData.h
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

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_H
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"
#include "MGIS/Behaviour/State.h"
#include "MGIS/Behaviour/Behaviour.h"
#ifdef __cplusplus
#include "MGIS/Behaviour/BehaviourData.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
/*! a simple alias */
using mgis_bv_BehaviourData = mgis::behaviour::BehaviourData;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_BehaviourData mgis_bv_BehaviourData;
#endif /*  __cplusplus */

/*!
 * \brief allocate behaviour data from a behaviour
 * \param[out] d: pointer to a pointer to the allocated data
 * \param[in]  b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_allocate_behaviour_data(
    mgis_bv_BehaviourData**, const mgis_bv_Behaviour* const);
/*!
 * \brief update the behaviour data
 * \param[out] d: pointer to a pointer to the allocated data
 */
MGIS_C_EXPORT mgis_status
mgis_bv_update_behaviour_data(mgis_bv_BehaviourData* const);
/*!
 * \brief reset the behaviour data
 * \param[out] d: pointer to a pointer to the allocated data
 */
MGIS_C_EXPORT mgis_status
mgis_bv_reset_behaviour_data(mgis_bv_BehaviourData* const);
/*!
 * \return the state at the beginning of the time step
 * \param[out] s: pointer to a pointer to the state
 * \param[in]  d: behaviour data
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_data_get_state_0(
    mgis_bv_State**, mgis_bv_BehaviourData* const);
/*!
 * \return the state at the end of the time step
 * \param[out] s: pointer to a pointer to the state
 * \param[in]  d: behaviour data
 */
MGIS_C_EXPORT mgis_status mgis_bv_behaviour_data_get_state_1(
    mgis_bv_State**, mgis_bv_BehaviourData* const);
/*!
 * \brief free the behaviour data
 * \param[out] d: pointer to a pointer to the allocated data
 */
MGIS_C_EXPORT mgis_status mgis_bv_free_behaviour_data(mgis_bv_BehaviourData**);


#ifdef __cplusplus
}
#endif /*  __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_H */
