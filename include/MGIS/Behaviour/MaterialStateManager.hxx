/*!
 * \file   MaterialStateManager.hxx
 * \brief
 * \author Thomas Helfer
 * \date   05/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX
#define LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX

#include <map>
#include <string>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Variant.hxx"
#include "MGIS/Span.hxx"
#include "MGIS/StringView.hxx"

namespace mgis {

  namespace behaviour {

    // forward declaration
    struct Behaviour;

    /*!
     * \brief structure in charge of handling the state of a
     * material in an optimized way. Here, the "material" is defined by a
     * behaviour and a number of integration points.
     *
     * The following design choices were made:
     * - The material properties and the external state variables are treated
     *   individually. They can be uniform or spatially variable.
     * - The internal state variables are treated as a block.
     */
    struct MGIS_EXPORT MaterialStateManager {
      //! a simple alias
      using FieldHolder =
          mgis::variant<real, mgis::span<real>, std::vector<real>>;
      /*!
       * \brief storage option for a non uniform material property or non
       * uniform external state variable.
       */
      enum StorageMode {
        LOCAL_STORAGE,    //! \brief use `std::vector` to store the data
        EXTERNAL_STORAGE  //! \brief use `mgis::span`  to store the data
      };                  // end of StorageMode
      /*!
       * \param[in] behaviour: behaviour
       * \param[in] s: number of integration points
       */
      MaterialStateManager(const Behaviour&, const size_type);
      //! move constructor
      MaterialStateManager(MaterialStateManager&&);
      //! copy constructor
      MaterialStateManager(const MaterialStateManager&);
      //! move assignement
      MaterialStateManager& operator=(MaterialStateManager&&);
      //! copy assignement
      MaterialStateManager& operator=(const MaterialStateManager&);
      //! \brief destructor
      ~MaterialStateManager();
      //! \brief value of the gradients
      std::vector<real> gradients;
      //! stride associate with the gradients
      const size_type gradients_stride;
      //! \brief values of the thermodynamic_forces
      std::vector<real> thermodynamic_forces;
      //! stride associate with the thermodynamic forces
      const size_type thermodynamic_forces_stride;
      /*!
       * \brief stored energy.
       */
      std::vector<real> stored_energies;
      /*!
       * \brief dissipated energies
       */
      std::vector<real> dissipated_energies;
      /*!
       * \brief material properties
       * The material properties can be uniform or not.
       * In the non uniform case, the data can be hold by the structure
       * (std::vector<real>) or simply borrow a reference (mgis::span<real>
       * case).
       */
      std::map<std::string, FieldHolder> material_properties;
      //! \brief values of the internal state variables
      std::vector<real> internal_state_variables;
      /*!
       * \brief stride associate with internal state variables.
       * \note this is also the size of an array containing all the internal
       * state variables for one integration point.
       */
      const size_type internal_state_variables_stride;
      /*!
       * \brief values of the external state variables
       * The external state variables can be uniform or not.
       * In the non uniform case, the data can be hold by the structure
       * (std::vector<real>) or simply borrow a reference (mgis::span<real>
       * case).
       */
      std::map<std::string, FieldHolder> external_state_variables;
      //! \brief number of integration points
      const size_type n;
      //! underlying behaviour
      const Behaviour& b;
    };  // end of struct MaterialStateManager

    /*!
     * \brief set the given material property
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: value
     */
    MGIS_EXPORT void setMaterialProperty(MaterialStateManager&,
                                         const mgis::string_view&,
                                         const real);
    /*!
     * \brief set the given material property
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT void setMaterialProperty(
        MaterialStateManager&,
        const mgis::string_view&,
        const mgis::span<real>&,
        const MaterialStateManager::StorageMode =
            MaterialStateManager::LOCAL_STORAGE);
    /*!
     * \return true if the given external state variable is defined.
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT bool isMaterialPropertyDefined(const MaterialStateManager&,
                                               const mgis::string_view&);
    /*!
     * \brief set the given material property
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT bool isMaterialPropertyUniform(const MaterialStateManager&,
                                               const mgis::string_view&);
    /*!
     * \return the uniform material property with the given name
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT real& getUniformMaterialProperty(MaterialStateManager&,
                                                 const mgis::string_view&);
    /*!
     * \return the uniform material property with the given name
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT const real& getUniformMaterialProperty(
        const MaterialStateManager&, const mgis::string_view&);
    /*!
     * \return the values of the material property with the given name
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT mgis::span<real> getNonUniformMaterialProperty(
        MaterialStateManager&, const mgis::string_view&);
    /*!
     * \return the values of the material property with the given name
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT mgis::span<const real> getNonUniformMaterialProperty(
        const MaterialStateManager&, const mgis::string_view&);
    /*!
     * \brief set the given external state variable
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: value
     */
    MGIS_EXPORT void setExternalStateVariable(MaterialStateManager&,
                                              const mgis::string_view&,
                                              const real);
    /*!
     * \brief set the given external state variable
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT void setExternalStateVariable(
        MaterialStateManager&,
        const mgis::string_view&,
        const mgis::span<real>&,
        const MaterialStateManager::StorageMode =
            MaterialStateManager::LOCAL_STORAGE);
    /*!
     * \return true if the given external state variable is defined.
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT bool isExternalStateVariableDefined(const MaterialStateManager&,
                                                    const mgis::string_view&);
    /*!
     * \return true if the given external state variable is uniform.
     * \param[out] m: material data manager
     * \param[in] n: name
     * \param[in] v: values
     * \param[in] s: storage mode
     */
    MGIS_EXPORT bool isExternalStateVariableUniform(const MaterialStateManager&,
                                                    const mgis::string_view&);
    /*!
     * \return the value of the uniform external state variable with the given name.
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT real& getUniformExternalStateVariable(MaterialStateManager&,
                                                      const mgis::string_view&);
    /*!
     * \return the uniform external state variable with the given name
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT const real& getUniformExternalStateVariable(
        const MaterialStateManager&, const mgis::string_view&);
    /*!
     * \return the values of the external state variable with the given name.
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT mgis::span<real> getNonUniformExternalStateVariable(
        MaterialStateManager&, const mgis::string_view&);
    /*!
     * \return the values of the external state variable with the given name.
     * \param[in] m: material data manager
     * \param[in] n: name
     */
    MGIS_EXPORT mgis::span<const real> getNonUniformExternalStateVariable(
        const MaterialStateManager&, const mgis::string_view&);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX */
