/*!
 * \file   Integrate.h
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

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_H
#define LIB_MGIS_BEHAVIOUR_INTEGRATE_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"
#include "MGIS/ThreadPool.h"
#include "MGIS/Behaviour/BehaviourDataView.h"
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/MaterialDataManager.h"

#ifdef __cplusplus
#include "MGIS/Behaviour/Integrate.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_bv_BehaviourIntegrationOptions =
    mgis::behaviour::BehaviourIntegrationOptions;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_bv_BehaviourIntegrationOptions
    mgis_bv_BehaviourIntegrationOptions;
#endif


//! kinematic of the behaviour treated
typedef enum {
  MGIS_BV_PREDICTION_TANGENT_OPERATOR = -3,
  MGIS_BV_PREDICTION_SECANT_OPERATOR = -2,
  MGIS_BV_PREDICTION_ELASTIC_OPERATOR = -1,
  MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR = 0,
  MGIS_BV_INTEGRATION_ELASTIC_OPERATOR = 1,
  MGIS_BV_INTEGRATION_SECANT_OPERATOR = 2,
  MGIS_BV_INTEGRATION_TANGENT_OPERATOR = 3,
  MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
} mgis_bv_IntegrationType;

//! integration options speed of sound flag
typedef enum {
  MGIS_BV_INTEGRATION_WITHOUT_SPEED_OF_SOUND = 0,
  MGIS_BV_INTEGRATION_WITH_SPEED_OF_SOUND = 1
} mgis_bv_SpeedOfSoundFlag;



/*!
 * \brief allocate a new
 * `mgis_bv_BehaviourIntegrationOptions` object
 *
 * \param[in] o: behaviour integration options
 */
MGIS_C_EXPORT mgis_status mgis_bv_create_behaviour_integration_options(
    mgis_bv_BehaviourIntegrationOptions**);
/*!
 * \brief set integration type
 * \param[in] o: behaviour integration options
 * \param[in] s: type
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_integration_options_set_integration_type(
    mgis_bv_BehaviourIntegrationOptions* const,
    const mgis_bv_IntegrationType);
/*!
 * \brief set speed of sound flag
 * \param[in] o: behaviour integration options
 * \param[in] s: flag
 */
MGIS_C_EXPORT mgis_status
mgis_bv_behaviour_integration_options_set_speed_of_sound_flag(
    mgis_bv_BehaviourIntegrationOptions* const,
    const mgis_bv_SpeedOfSoundFlag);
/*!
 * \brief free the ressources associated with a
 * `mgis_bv_BehaviourIntegrationOptions` object
 *
 * \param[in] o: integration behaviour options
 */
MGIS_C_EXPORT mgis_status mgis_bv_free_behaviour_integration_options(
    mgis_bv_BehaviourIntegrationOptions**);
/*!
 * \brief integrate the behaviour. The returned value has the following
 * meaning:
 * - -1: integration failed
 * -  0: integration succeeded but results are unreliable
 * -  1: integration succeeded and results are reliable
 *
 * \param[out] r: result
 * \param[in,out] d: behaviour data
 * \param[in,out] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_integrate(int* const,
                                            mgis_bv_BehaviourDataView* const,
                                            const mgis_bv_Behaviour* const);
//! \brief debugging version
MGIS_C_EXPORT mgis_status
mgis_bv_integrate_debug(int* const,
                        mgis_bv_BehaviourDataView* const,
                        const mgis_bv_Behaviour* const);
/*!
 * \brief integrate the behaviour. The returned value has the following
 * meaning:
 * - -1: integration failed
 * -  0: integration succeeded but results are unreliable
 * -  1: integration succeeded and results are reliable
 *
 * \note this version internally creates a view of the data and calls
 * `mgis_bv_integrate`
 * \note this function has been introduced mostly for bindings to
 * languages where exposing `mgis_bv_BehaviourDataView` would be too
 * difficult (see the `fortran` bindings for example).
 * \param[out] r: result
 * \param[in,out] d: behaviour data
 * \param[in,out] b: behaviour
 */
MGIS_C_EXPORT mgis_status mgis_bv_integrate_2(int* const,
					      mgis_bv_BehaviourData* const,
					      const mgis_bv_Behaviour* const);
//! \brief debugging version
MGIS_C_EXPORT mgis_status mgis_bv_integrate_debug_2(
    int* const, mgis_bv_BehaviourData* const, const mgis_bv_Behaviour* const);
/*!
 * \brief integrate the behaviour for a range of integration points. The
 * returned value has the following meaning:
 * - -1: integration failed for at least one Gauss point
 * -  0: all integrations succeeded but results are unreliable for at least
 *       one Gauss point
 * -  1: integration succeeded and results are reliable.
 *
 * \param[out] r: result
 * \param[in,out] m: material data manager
 * \param[in] i: type of integration to be performed
 * \param[in] dt: time step
 * \param[in] b: first index of the range
 * \param[in] e: last index of the range
 */
MGIS_C_EXPORT mgis_status
mgis_bv_integrate_material_data_manager_part(int* const,
                                             mgis_bv_MaterialDataManager* const,
                                             const mgis_bv_IntegrationType,
                                             const mgis_real,
                                             const mgis_size_type,
                                             const mgis_size_type);
//! \brief debugging version
MGIS_C_EXPORT mgis_status
mgis_bv_integrate_debug_material_data_manager_part(int* const,
                                             mgis_bv_MaterialDataManager* const,
                                             const mgis_bv_IntegrationType,
                                             const mgis_real,
                                             const mgis_size_type,
                                             const mgis_size_type);
/*!
 * \brief integrate the behaviour for a range of integration points. The
 * returned value has the following meaning:
 * - -1: integration failed for at least one Gauss point
 * -  0: all integrations succeeded but results are unreliable for at least
 *       one Gauss point
 * -  1: integration succeeded and results are reliable.
 *
 * \param[out] r: result
 * \param[in,out] p: thread pool
 * \param[in,out] m: material data manager
 * \param[in] i: type of integration to be performed
 * \param[in] dt: time step
 */
MGIS_C_EXPORT mgis_status
mgis_bv_integrate_material_data_manager(int* const,
                                        mgis_ThreadPool* const,
                                        mgis_bv_MaterialDataManager* const,
                                        const mgis_bv_IntegrationType,
                                        const mgis_real);
//! \brief debugging version
MGIS_C_EXPORT mgis_status mgis_bv_integrate_debug_material_data_manager(
    int* const,
    mgis_ThreadPool* const,
    mgis_bv_MaterialDataManager* const,
    const mgis_bv_IntegrationType,
    const mgis_real);

MGIS_C_EXPORT mgis_status
mgis_bv_integrate_material_data_manager_with_options(int* const,
                                                     mgis_ThreadPool* const,
                                                     mgis_bv_MaterialDataManager* const,
                                                     mgis_bv_BehaviourIntegrationOptions* const,
                                                     const mgis_real);
//! \brief debugging version
MGIS_C_EXPORT mgis_status
mgis_bv_integrate_debug_material_data_manager_with_options(
    int* const,
    mgis_ThreadPool* const,
    mgis_bv_MaterialDataManager* const,
    mgis_bv_BehaviourIntegrationOptions* const,
    const mgis_real);

#ifdef __cplusplus
}
#endif /*  __cplusplus */

#endif /*LIB_MGIS_BEHAVIOUR_INTEGRATE_H */
