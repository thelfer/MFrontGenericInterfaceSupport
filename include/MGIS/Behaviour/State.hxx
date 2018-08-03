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
      //! default constructor
      State();
      //! move constructor
      State(State&&);
      //! copy constructor
      State(const State&);
      //! move assignement
      State& operator=(State&&);
      //! copy assignement
      State& operator=(const State&);
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
     * \brief set the value of a material property
     * \param[out] d: behaviour data
     * \param[in]  b: behaviour
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setMaterialProperty(State&,
                                         const Behaviour&,
                                         const string_view,
                                         const real);
    /*!
     * \brief set the value of a material property
     * \param[out] d: behaviour data
     * \param[in]  o: material property offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setMaterialProperty(State&,
                                         const size_type,
                                         const real);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] d: behaviour data
     * \param[in]  b: behaviour
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const Behaviour&,
                                              const string_view,
                                              const real);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] d: behaviour data
     * \param[in]  b: behaviour
     * \param[in]  n: name
     * \param[in]  v: values
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const Behaviour&,
                                              const string_view,
                                              const real* const);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] d: behaviour data
     * \param[in]  o: internal state variable offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const size_type,
                                              const real);
    /*!
     * \brief set the value of an internal state variable
     * \param[out] d: behaviour data
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
     * \param[out] d: behaviour data
     * \param[in]  o: internal state variable offset
     * \param[in]  n: internal state variable size
     * \param[in]  v: values
     */
    MGIS_EXPORT void setInternalStateVariable(State&,
                                              const size_type,
                                              const size_type,
                                              const real* const);
    /*!
     * \brief set the value of a external state variable
     * \param[out] d: behaviour data
     * \param[in]  b: behaviour
     * \param[in]  n: name
     * \param[in]  v: value
     */
    MGIS_EXPORT void setExternalStateVariable(State&,
                                              const Behaviour&,
                                              const string_view,
                                              const real);
    /*!
     * \brief set the value of a external state variable
     * \param[out] d: behaviour data
     * \param[in]  o: external state variable offset
     * \param[in]  v: value
     */
    MGIS_EXPORT void setExternalStateVariable(State&,
                                              const size_type,
                                              const real);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_STATE_HXX */
