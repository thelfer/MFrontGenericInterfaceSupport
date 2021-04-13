/*!
 * \file   include/MGIS/Behaviour/BehaviourDataView.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/07/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_HXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_HXX

#ifdef __cplusplus
#include "MGIS/Config.hxx"
#else
#include "MGIS/Config-c.h"
#endif /* __cplusplus */

#include "MGIS/Behaviour/StateView.hxx"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct mgis_bv_BehaviourDataView mgis_bv_BehaviourDataView;

/*!
 * \brief structure
 */
struct mgis_bv_BehaviourDataView {
  /*!
   * \brief time increment
   */
  mgis_real dt;
  /*!
   * \brief the stiffness matrix.
   * On input, the first element of K must contain the type of type
   * of stiffness matrix expected. If this value is negative, only
   * the prediction operator is computed. This value has the
   * following meaning:
   * - if K[0] is lower than -2.5, the tangent operator must be
   *   computed.
   * - if K[0] is in [-2.5:-1.5]: the secant operator must be
   *   computed.
   * - if K[0] is in [-1.5:-0.5]: the elastic operator must be
   *   computed.
   * - if K[0] is in [-0.5:0.5]: the behaviour integration is
   *   performed, but no stiffness matrix.
   * - if K[0] is in [0.5:1.5]: the elastic operator must be
   *   computed.
   * - if K[0] is in [1.5:2.5]: the secant operator must be
   *   computed.
   * - if K[0] is in [2.5:3.5]: the secant operator must be
   *   computed.
   * - if K[0] is greater than 3.5, the consistent tangent operator
   *   must be computed.
   *
   * Other values of K are meant to store behaviour's option. This
   * is currently only meaningful for finite strain behaviours.
   *
   * For finite strain behaviours, K[1] holds the choice of the stress
   * measure used:
   * - if K[1] < 0.5, the Cauchy stress is used
   * - if 0.5 < K[1] < 1.5, the second Piola-Kirchoff stress is used
   * - if 1.5 < K[1] < 2.5, the first Piola-Kirchoff stress is used
   *
   * For finite strain behaviours, K[2] holds the choice of the
   * consistent tangent operator returned by the behaviour:
   * - if K[2]<0.5, the derivative of the Cauchy stress with respect
   *   to the deformation gradient is returned
   * - if 0.5<K[2]<1.5, the derivative of the second Piola-Kirchoff
   *   stress with respect to the Green-Lagrange strain
   *   is returned
   * - if 1.5<K[2]<2.5, the derivative of the first Piola-Kirchoff
   *   stress with respect to the deformation gradient is returned
   */
  mgis_real* K;
  /*!
   * \brief proposed time step increment increase factor
   */
  mgis_real rdt;
  /*!
   * \brief state at the beginning of the time step
   */
  mgis_bv_InitialStateView s0;
  /*!
   * \brief state at the end of the time step
   */
  mgis_bv_StateView s1;
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus

namespace mgis {

  namespace behaviour {

    //! a simple alias
    using BehaviourDataView = ::mgis_bv_BehaviourDataView;

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATAVIEW_HXX */
