/*!
 * \file   State.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_STATE_HXX
#define LIB_MGIS_BEHAVIOUR_STATE_HXX

#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/StringView.hxx"
#include "MGIS/Behaviour/StateView.hxx"

namespace mgis {

  namespace behaviour {

    // forward declaration
    struct Behaviour;

    /*!
     * \brief structure in charge of containing the state of a material.
     */
    struct MGIS_EXPORT State {
      /*!
       * \brief constructor from a behaviour
       * \param[in] b: behaviour
       */
      State(const Behaviour&);
      //! move constructor
      State(State&&);
      //! copy constructor
      State(const State&);
      //! move assignement
      State& operator=(State&&);
      //! copy assignement
      State& operator=(const State&);
      //! \brief the underlying behaviour
      const Behaviour& b;
      /*!
       * \brief stored energy (computed by `@InternalEnergy` in `MFront`
       * files) This output is optional.
       */
      real stored_energy;
      /*!
       * \brief dissipated energy (computed by `@DissipatedEnergy` in `MFront`
       * files) This output is optional.
       */
      real dissipated_energy;
      //! \brief value of the gradients
      std::vector<real> gradients;
      //! \brief values of the thermodynamic_forces
      std::vector<real> thermodynamic_forces;
      //! \brief values of the material properties
      std::vector<real> material_properties;
      //! \brief values of the internal state variables
      std::vector<real> internal_state_variables;
      //! \brief values of the external state variables
      std::vector<real> external_state_variables;
    };  // end of struct State

    /*!
     * \brief set the value of a gradient
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setGradient(State&, const string_view, const real);
    /*!
     * \brief set the value of a gradient
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: values
     */
    MGIS_EXPORT void setGradient(State&, const string_view, const real* const);
    /*!
     * \brief set the value of a gradient
     * \param[out] s: state
     * \param[in]  o: gradient offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setGradient(State&, const size_type, const real);
    /*!
     * \brief set the value of a gradient
     * \param[out] s: state
     * \param[in]  o: gradient offset
     * \param[in]  n: gradient size
     * \param[in]  v: value
     */
    MGIS_EXPORT void setGradient(State&,
                                 const size_type,
                                 const size_type,
                                 const real);
    /*!
     * \brief set the values of a gradient
     * \param[out] s: state
     * \param[in]  o: gradient offset
     * \param[in]  n: gradient size
     * \param[in]  v: values
     */
    MGIS_EXPORT void setGradient(State&,
                                 const size_type,
                                 const size_type,
                                 const real* const);
    /*!
     * \brief get the value(s) of a gradient
     * \return a pointer the value(s) of the gradient
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT real* getGradient(State&, const string_view);
    /*!
     * \brief get the value(s) of a gradient
     * \return a pointer the value(s) of the gradient
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT const real* getGradient(const State&, const string_view);
    /*!
     * \brief get the value(s) of a gradient
     * \return a pointer the value(s) of the gradient
     * \param[in] s: state
     * \param[in] o: gradient offset
     */
    MGIS_EXPORT real* getGradient(State&, const size_type);
    /*!
     * \brief get the value(s) of a gradient
     * \return a pointer the value(s) of the gradient
     * \param[in] s: state
     * \param[in] o: gradient offset
     */
    MGIS_EXPORT const real* getGradient(const State&, const size_type);
    /*!
     * \brief set the value of a thermodynamic force
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setThermodynamicForce(State&,
                                           const string_view,
                                           const real);
    /*!
     * \brief set the value of a thermodynamic force
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: values
     */
    MGIS_EXPORT void setThermodynamicForce(State&,
                                           const string_view,
                                           const real* const);
    /*!
     * \brief set the value of a thermodynamic force
     * \param[out] s: state
     * \param[in]  o: thermodynamic force offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setThermodynamicForce(State&, const size_type, const real);
    /*!
     * \brief set the value of a thermodynamic force
     * \param[out] s: state
     * \param[in]  o: thermodynamic force offset
     * \param[in]  n: thermodynamic force size
     * \param[in]  v: value
     */
    MGIS_EXPORT void setThermodynamicForce(State&,
                                           const size_type,
                                           const size_type,
                                           const real);
    /*!
     * \brief set the values of a thermodynamic force
     * \param[out] s: state
     * \param[in]  o: thermodynamic force offset
     * \param[in]  n: thermodynamic force size
     * \param[in]  v: values
     */
    MGIS_EXPORT void setThermodynamicForce(State&,
                                           const size_type,
                                           const size_type,
                                           const real* const);
    /*!
     * \brief get the value(s) of a thermodynamic force
     * \return a pointer the value(s) of the thermodynamic force
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT real* getThermodynamicForce(State&, const string_view);
    /*!
     * \brief get the value(s) of a thermodynamic force
     * \return a pointer the value(s) of the thermodynamic force
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT const real* getThermodynamicForce(const State&,
                                                  const string_view);
    /*!
     * \brief get the value(s) of a thermodynamic force
     * \return a pointer the value(s) of the thermodynamic force
     * \param[in] s: state
     * \param[in] o: thermodynamic force offset
     */
    MGIS_EXPORT real* getThermodynamicForce(State&, const size_type);
    /*!
     * \brief get the value(s) of a thermodynamic force
     * \return a pointer the value(s) of the thermodynamic force
     * \param[in] s: state
     * \param[in] o: thermodynamic force offset
     */
    MGIS_EXPORT const real* getThermodynamicForce(const State&,
                                                  const size_type);
    /*!
     * \brief set the value of a material property
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setMaterialProperty(State&, const string_view, const real);
    /*!
     * \return a pointer to the value of a material property
     * \param[in] s: state
     * \param[in] n: name
     */
    MGIS_EXPORT real* getMaterialProperty(State&, const string_view);
    /*!
     * \return a pointer to the value of a material property
     * \param[in] s: state
     * \param[in] n: name
     */
    MGIS_EXPORT const real* getMaterialProperty(const State&,
                                                const string_view);
    /*!
     * \brief set the value of a material property
     * \param[out] s: state
     * \param[in]  o: material property offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setMaterialProperty(State&, const size_type, const real);
    /*!
     * \return a pointer to the value of a material property
     * \param[out] s: state
     * \param[in]  o: material property offset
     */
    MGIS_EXPORT real* getMaterialProperty(State&, const size_type);
    /*!
     * \return a pointer to the value of a material property
     * \param[out] s: state
     * \param[in]  o: material property offset
     */
    MGIS_EXPORT const real* getMaterialProperty(const State&, const size_type);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const string_view,
                                              const real);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: values
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const string_view,
                                              const real* const);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] s: state
     * \param[in]  o: internal state variable offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const size_type,
                                              const real);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] s: state
     * \param[in]  o: internal state variable offset
     * \param[in]  n: internal state variable size
     * \param[in]  v: value
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const size_type,
                                              const size_type,
                                              const real);
    /*!
     * \brief set the values of an internal state variable
     * \param[out] s: state
     * \param[in]  o: internal state variable offset
     * \param[in]  n: internal state variable size
     * \param[in]  v: values
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const size_type,
                                              const size_type,
                                              const real* const);
    /*!
     * \brief get the value(s) of an internal state variable
     * \return a pointer the value(s) of the internal state variable
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT real* getInternalStateVariable(State&, const string_view);
    /*!
     * \brief get the value(s) of an internal state variable
     * \return a pointer the value(s) of the internal state variable
     * \param[in] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT const real* getInternalStateVariable(const State&,
                                                     const string_view);
    /*!
     * \brief get the value(s) of an internal state variable
     * \return a pointer the value(s) of the internal state variable
     * \param[in] s: state
     * \param[in] o: internal state variable offset
     */
    MGIS_EXPORT real* getInternalStateVariable(State&, const size_type);
    /*!
     * \brief get the value(s) of an internal state variable
     * \return a pointer the value(s) of the internal state variable
     * \param[in] s: state
     * \param[in] o: internal state variable offset
     */
    MGIS_EXPORT const real* getInternalStateVariable(const State&,
                                                     const size_type);
    /*!
     * \brief set the value of an external state variable
     * \param[out] s: state
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setExternalStateVariable(State&,
                                              const string_view,
                                              const real);
    /*!
     * \brief set the value of an external state variable
     * \param[out] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT real* getExternalStateVariable(State&, const string_view);
    /*!
     * \brief set the value of an external state variable
     * \param[out] s: state
     * \param[in]  n: name
     */
    MGIS_EXPORT const real* getExternalStateVariable(const State&,
                                                     const string_view);
    /*!
     * \brief set the value of an external state variable
     * \param[out] s: state
     * \param[in]  o: external state variable offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setExternalStateVariable(State&,
                                              const size_type,
                                              const real);
    /*!
     * \return a pointer to the value of an external state variable
     * \param[out] s: state
     * \param[in]  o: external state variable offset
     */
    MGIS_EXPORT real* getExternalStateVariable(State&, const size_type);
    /*!
     * \return a pointer to the value of an external state variable
     * \param[out] s: state
     * \param[in]  o: external state variable offset
     */
    MGIS_EXPORT const real* getExternalStateVariable(const State&,
                                                     const size_type);

    /*!
     * \brief make a view from a behaviour data
     * \param[in] s: state
     * \return the view
     * \note the view has at most the same life time as the state.
     */
    MGIS_EXPORT StateView make_view(State&);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_STATE_HXX */
