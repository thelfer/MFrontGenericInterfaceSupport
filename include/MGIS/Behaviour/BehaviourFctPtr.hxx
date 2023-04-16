/*!
 * \file   BehaviourFctPtr.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX

#include "MGIS/Config-c.h"

// forward declaration
typedef struct mgis_bv_BehaviourDataView mgis_bv_BehaviourDataView;
//! \brief type of a pointer function implementing the behaviour integration
typedef int (*mgis_bv_BehaviourFctPtr)(mgis_bv_BehaviourDataView* const);
/*!
 * \brief type of the pointer of a function implementing the rotation of the
 * gradients from the global frame to the material frame
 * \param[in] mg: gradients in the material frame
 * \param[in] gg: gradients in the global frame
 * \param[in] r: rotation matrix
 */
typedef void (*mgis_bv_RotateBehaviourGradientsFctPtr)(mgis_real* const,
                                                       const mgis_real* const,
                                                       const mgis_real* const);
/*!
 * \brief type of the pointer of a function implementing the rotation of the
 * gradients from the global frame to the material frame
 */
typedef void (*mgis_bv_RotateArrayOfBehaviourGradientsFctPtr)(
    mgis_real* const,
    const mgis_real* const,
    const mgis_real* const,
    const mgis_size_type);
/*!
 * \brief type of the pointer of a function implementing the rotation of the
 * thermodynamic forces from the global frame to the material frame
 */
typedef void (*mgis_bv_RotateBehaviourThermodynamicForcesFctPtr)(
    mgis_real* const, const mgis_real* const, const mgis_real* const);
/*!
 * \brief type of the pointer of a function implementing the rotation of an
 * array of thermodynamic forces from the global frame to the material frame
 */
typedef void (*mgis_bv_RotateArrayOfBehaviourThermodynamicForcesFctPtr)(
    mgis_real* const,
    const mgis_real* const,
    const mgis_real* const,
    const mgis_size_type);
/*!
 * \brief type of the pointer of a function implementing the rotation of the
 * tangent operator blocks from the global frame to the material frame
 */
typedef void (*mgis_bv_RotateBehaviourTangentOperatorBlocksFctPtr)(
    mgis_real* const, const mgis_real* const, const mgis_real* const);
/*!
 * \brief type of the pointer of a function implementing the rotation of an
 * array of tangent operator blocks from the global frame to the material frame
 */
typedef void (*mgis_bv_RotateArrayOfBehaviourTangentOperatorBlocksFctPtr)(
    mgis_real* const,
    const mgis_real* const,
    const mgis_real* const,
    const mgis_size_type);

#ifdef __cplusplus

namespace mgis::behaviour {

  //! \brief a simple alias
  using BehaviourFctPtr = mgis_bv_BehaviourFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of the
   * gradients from the global frame to the material frame
   */
  using RotateBehaviourGradientsFctPtr = mgis_bv_RotateBehaviourGradientsFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of an
   * array of gradients from the global frame to the material frame
   */
  using RotateArrayOfBehaviourGradientsFctPtr =
      mgis_bv_RotateArrayOfBehaviourGradientsFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of the
   * thermodynamic forces from the global frame to the material frame
   */
  using RotateBehaviourThermodynamicForcesFctPtr =
      mgis_bv_RotateBehaviourThermodynamicForcesFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of an
   * array of thermodynamic forces from the global frame to the material frame
   */
  using RotateArrayOfBehaviourThermodynamicForcesFctPtr =
      mgis_bv_RotateArrayOfBehaviourThermodynamicForcesFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of the
   * tangent operator blocks from the global frame to the material frame.
   */
  using RotateBehaviourTangentOperatorBlocksFctPtr =
      mgis_bv_RotateBehaviourTangentOperatorBlocksFctPtr;
  /*!
   * \brief type of the pointer of a function implementing the rotation of an
   * array of tangent operator blocks from the global frame to the material
   * frame.
   */
  using RotateArrayOfBehaviourTangentOperatorBlocksFctPtr =
      mgis_bv_RotateArrayOfBehaviourTangentOperatorBlocksFctPtr;

}  // end of namespace mgis::behaviour

#endif /* __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURFCTPTR_HXX */
