/*!
 * \file   include/MGIS/Behaviour/MaterialStateManager.hxx
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
#include <variant>
#include <optional>
#include <string_view>
#include "MGIS/Config.hxx"
#include "MGIS/Span.hxx"
#include "MGIS/StorageMode.hxx"

namespace mgis::behaviour {

  // forward declaration
  struct Behaviour;

  /*!
   * \brief a structure in charge of holding information on how a material
   * state manager shall be initialized.
   * It may contain pointers to externally allocated data, that won't be
   * handled by the final state manager.
   * If a pointer is not initialized, the material state manager will allocate
   * and handle memory internally.
   */
  struct MaterialStateManagerInitializer {
    /*!
     * \brief view to an externally allocated memory used to store the
     * gradients. If empty, the material state manager will
     * initialize the required memory internally.
     */
    mgis::span<mgis::real> gradients;
    /*!
     * \brief view to an externally allocated memory used to store the
     * thermodynamic forces. If empty, the material state manager will
     * initialize the required memory internally.
     */
    mgis::span<mgis::real> thermodynamic_forces;
    /*!
     * \brief view to an externally allocated memory used to store the
     * internal state variables. If empty, the material state manager will
     * initialize the required memory internally.
     */
    mgis::span<mgis::real> internal_state_variables;
    /*!
     * \brief view to an externally allocated memory used to store the stored
     * energies. If empty, the material state manager will initialize the
     * required memory internally, if the behaviour computes the stored
     * energy.
     *
     * \note for backward compatibililty, the user may allocate memory for the
     * stored energies even if the behaviour don't compute the stored
     * energy.
     */
    mgis::span<mgis::real> stored_energies;
    /*!
     * \brief view to an externally allocated memory used to store the
     * dissipated energies. If empty, the material state manager will
     * initialize the required memory internally, if the behaviour computes
     * the dissipated energy.
     *
     * \note for backward compatibililty, the user may allocate memory for the
     * dissipated energies even if the behaviour don't compute the dissipated
     * energy.
     */
    mgis::span<mgis::real> dissipated_energies;
  };  // end of MaterialStateManagerInitializer

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
    //! \brief a simple alias
    using FieldHolder =
        std::variant<real, mgis::span<mgis::real>, std::vector<mgis::real>>;
    //! \brief a simple alias
    using StorageMode = mgis::StorageMode;
    //!
    static constexpr auto LOCAL_STORAGE = mgis::StorageMode::LOCAL_STORAGE;
    //!
    static constexpr auto EXTERNAL_STORAGE =
        mgis::StorageMode::EXTERNAL_STORAGE;
    /*!
     * \param[in] behaviour: behaviour
     * \param[in] s: number of integration points
     */
    MaterialStateManager(const Behaviour&, const size_type);
    /*!
     * \param[in] behaviour: behaviour
     * \param[in] s: number of integration points
     * \param[in] i: initializer
     */
    MaterialStateManager(const Behaviour&,
                         const size_type,
                         const MaterialStateManagerInitializer&);
    //! \brief destructor
    ~MaterialStateManager();
    //! \brief view to the values of the gradients
    mgis::span<mgis::real> gradients;
    //! stride associate with the gradients
    const size_type gradients_stride;
    //! \brief view to the values of the thermodynamic_forces
    mgis::span<mgis::real> thermodynamic_forces;
    //! stride associate with the thermodynamic forces
    const size_type thermodynamic_forces_stride;
    /*!
     * \brief view to the values to the stored energy.
     */
    mgis::span<mgis::real> stored_energies;
    /*!
     * \brief view to the values to the dissipated energies
     */
    mgis::span<mgis::real> dissipated_energies;
    /*!
     * \brief material properties
     * The material properties can be uniform or not.
     * In the non uniform case, the data can be hold by the structure
     * (std::vector<real>) or simply borrow a reference
     * (mgis::span<mgis::real>
     * case).
     */
    std::map<std::string, FieldHolder> material_properties;
    /*! \brief mass density
     * The mass density can be uniform or not.
     * In the non uniform case, the data can be hold by the structure
     * (std::vector<real>) or simply borrow a reference
     * (mgis::span<mgis::real>
     * case).
     */
    std::optional<FieldHolder> mass_density;
    //! \brief view to the values of the internal state variables
    mgis::span<mgis::real> internal_state_variables;
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
     * (std::vector<real>) or simply borrow a reference
     * (mgis::span<mgis::real>
     * case).
     */
    std::map<std::string, FieldHolder> external_state_variables;
    //! \brief number of integration points
    const size_type n;
    //! underlying behaviour
    const Behaviour& b;

   private:
    //! \brief value of the gradients, if hold internally
    std::vector<mgis::real> gradients_values;
    //! \brief value of the thermodynamic forces, if hold internally
    std::vector<mgis::real> thermodynamic_forces_values;
    //! \brief value of the internal state variables, if hold internally
    std::vector<mgis::real> internal_state_variables_values;
    //! \brief value of the stored energies, if hold internally
    std::vector<mgis::real> stored_energies_values;
    //! \brief value of the dissipated energies, if hold internally
    std::vector<mgis::real> dissipated_energies_values;
    //! \brief move constructor
    MaterialStateManager(MaterialStateManager&&) = delete;
    //! \brief copy constructor
    MaterialStateManager(const MaterialStateManager&) = delete;
    //! \brief move assignement
    MaterialStateManager& operator=(MaterialStateManager&&) = delete;
    //! \brief copy assignement
    MaterialStateManager& operator=(const MaterialStateManager&) = delete;
  };  // end of struct MaterialStateManager

  /*!
   * \brief set the given material property
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: value
   */
  MGIS_EXPORT void setMaterialProperty(MaterialStateManager&,
                                       const std::string_view&,
                                       const real);
  /*!
   * \brief set the given material property
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: values
   * \param[in] s: storage mode
   */
  MGIS_EXPORT void setMaterialProperty(MaterialStateManager&,
                                       const std::string_view&,
                                       const mgis::span<mgis::real>&,
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
                                             const std::string_view&);
  /*!
   * \brief chek if the given material property is uniform
   * \param[out] m: material data manager
   * \param[in] n: name
   */
  MGIS_EXPORT bool isMaterialPropertyUniform(const MaterialStateManager&,
                                             const std::string_view&);
  /*!
   * \brief set the mass density
   * \param[out] m: material data manager
   * \param[in] v: value
   */
  MGIS_EXPORT void setMassDensity(MaterialStateManager&, const real);
  /*!
   * \brief set the mass density
   * \param[out] m: material data manager
   * \param[in] v: values
   * \param[in] s: storage mode
   */
  MGIS_EXPORT void setMassDensity(MaterialStateManager&,
                                  const mgis::span<mgis::real>&,
                                  const MaterialStateManager::StorageMode =
                                      MaterialStateManager::LOCAL_STORAGE);
  /*!
   * \return true if the given external state variable is defined.
   * \param[out] m: material data manager
   */
  MGIS_EXPORT bool isMassDensityDefined(const MaterialStateManager&);
  /*!
   * \return true if the mass density is uniform
   * \param[out] m: material data manager
   */
  MGIS_EXPORT bool isMassDensityUniform(const MaterialStateManager&);
  /*!
   * \brief set the given external state variable
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: value
   */
  MGIS_EXPORT void setExternalStateVariable(MaterialStateManager&,
                                            const std::string_view&,
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
      const std::string_view&,
      const mgis::span<mgis::real>&,
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
                                                  const std::string_view&);
  /*!
   * \return true if the given external state variable is uniform.
   * \param[out] m: material data manager
   * \param[in] n: name
   */
  MGIS_EXPORT bool isExternalStateVariableUniform(const MaterialStateManager&,
                                                  const std::string_view&);
  /*!
   * \brief update the values of a state from another state
   * \param[out] o: output state
   * \param[out] i: input state
   */
  MGIS_EXPORT void updateValues(MaterialStateManager&,
                                const MaterialStateManager&);
  /*!
   * \brief extract an internal state variable
   *
   * \param[out] o: buffer in which the values of the given internal state
   * variable is stored
   * \param[in] s: material state manager
   * \param[in] n: name of the internal state variables
   *
   * \note the output buffer must be allocated properly
   */
  MGIS_EXPORT void extractInternalStateVariable(
      mgis::span<mgis::real>,
      const mgis::behaviour::MaterialStateManager&,
      const std::string_view);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX */
