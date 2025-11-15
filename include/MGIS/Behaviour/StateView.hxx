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
#include <iosfwd>
#include "MGIS/Config.hxx"
#else
#include "MGIS/Config-c.h"
#endif /* __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

//! \brief state of the material at the end of the time step
typedef struct {
  //! \brief value of the gradients
  const mgis_real* gradients;
  //! \brief values of the thermodynamic_forces
  mgis_real* thermodynamic_forces;
  //! \brief volumetric mass density in the reference configuration
  const mgis_real* mass_density;
  //! \brief values of the material properties
  const mgis_real* material_properties;
  //! \brief values of the internal state variables
  mgis_real* internal_state_variables;
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
  const mgis_real* external_state_variables;
} mgis_bv_StateView;

//! \brief state of the material at the beginning of the time step
typedef struct {
  //! \brief value of the gradients
  const mgis_real* gradients;
  //! \brief values of the thermodynamic_forces
  const mgis_real* thermodynamic_forces;
  //! \brief volumetric mass density in the reference configuration
  const mgis_real* mass_density;
  //! \brief values of the material properties
  const mgis_real* material_properties;
  //! \brief values of the internal state variables
  const mgis_real* internal_state_variables;
  /*!
   * \brief stored energy (computed by `@InternalEnergy` in `MFront`
   * files) This output is optional.
   */
  const mgis_real* stored_energy;
  /*!
   * \brief stored energy (computed by `@DissipatedEnergy` in `MFront`
   * files) This output is optional.
   */
  const mgis_real* dissipated_energy;
  //! \brief values of the external state variables
  const mgis_real* external_state_variables;
} mgis_bv_InitialStateView;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus

namespace mgis::behaviour {

  // forward declaration
  struct Behaviour;

  //! a simple alias
  using StateView = ::mgis_bv_StateView;
  //! a simple alias
  using InitialStateView = ::mgis_bv_InitialStateView;

  /*!
   * \brief print a detailled (verbose) description of the integration point
   * state using a markdown format
   * \param[in] os: ouptut stream
   * \param[in] b: behaviour
   * \param[in] s: state
   * \param[in] l: title level
   */
  MGIS_EXPORT void print_markdown(std::ostream&,
                                  const Behaviour&,
                                  const StateView&,
                                  const mgis::size_type);
  /*!
   * \brief print a detailled (verbose) description of the integration point
   * state using a markdown format
   * \param[in] os: ouptut stream
   * \param[in] b: behaviour
   * \param[in] s: state
   * \param[in] l: title level
   */
  MGIS_EXPORT void print_markdown(std::ostream&,
                                  const Behaviour&,
                                  const InitialStateView&,
                                  const mgis::size_type);

}  // end of namespace mgis::behaviour

#endif /* __cplusplus */

#endif /* LIB_MGIS_BEHAVIOUR_STATEVIEW_HXX */
