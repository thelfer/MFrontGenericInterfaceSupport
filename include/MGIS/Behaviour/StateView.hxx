/*!
 * \file   include/MGIS/Behaviour/StateView.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_STATEVIEW_HXX
#define LIB_MGIS_BEHAVIOUR_STATEVIEW_HXX

#ifdef __cplusplus
#include "MGIS/Config.hxx"
#else
#include "MGIS/Config-c.h"
#endif /* __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
 * \brief state of the material
 */
typedef struct{
  //! \brief value of the gradients
  mgis_real * gradients;
  //! \brief values of the thermodynamic_forces
  mgis_real * thermodynamic_forces;
  //! \brief values of the material properties
  mgis_real * material_properties;
  //! \brief values of the internal state variables
  mgis_real * internal_state_variables;
  /*!
   * \brief stored energy (computed by `@InternalEnergy` in `MFront`
   * files) This output is optional.
   */
  mgis_real* stored_energy;
  /*!
   * \brief stored energy (computed by `@DissipatedEnergy` in `MFront`
   * files) This output is optional.
   */
  mgis_real* dissipated_energy;
  //! \brief values of the external state variables
  mgis_real * external_state_variables;
} mgis_bv_StateView;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus

namespace mgis {

  namespace behaviour {

    //! a simple alias
    using StateView = ::mgis_bv_StateView;

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_STATEVIEW_HXX */
