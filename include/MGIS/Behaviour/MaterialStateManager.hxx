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
#include <span>
#include <string>
#include <vector>
#include <variant>
#include <optional>
#include <string_view>
#ifdef MGIS_HAVE_HDF5
#include "MGIS/Utilities/HDF5Forward.hxx"
#endif /* MGIS_HAVE_HDF5 */

#include "MGIS/Config.hxx"
#include "MGIS/StorageMode.hxx"

#ifdef MGIS_FUNCTION_SUPPORT
#include "MGIS/Function/Function.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

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
    std::span<mgis::real> gradients;
    /*!
     * \brief view to an externally allocated memory used to store the
     * thermodynamic forces. If empty, the material state manager will
     * initialize the required memory internally.
     */
    std::span<mgis::real> thermodynamic_forces;
    /*!
     * \brief view to an externally allocated memory used to store the
     * internal state variables. If empty, the material state manager will
     * initialize the required memory internally.
     */
    std::span<mgis::real> internal_state_variables;
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
    std::span<mgis::real> stored_energies;
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
    std::span<mgis::real> dissipated_energies;
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
    struct FieldHolder {
      FieldHolder& operator=(const mgis::real) noexcept;
      //! \brief pointer to the values of the field
      std::variant<real, std::span<mgis::real>, std::vector<mgis::real>> value;
      /*!
       * \brief boolean stating if the values shall be updated by the
       * `updateValues` which is call by the revert or update functions.
       *
       * Setting this flag is usefull when the update of the field is already
       * handled.
       */
      bool shall_be_updated = true;
    };
    /*!
     * \brief enum used to express if a variable (material property, external
     * state variables, mass density) shall be updated or reverted automatically
     */
    enum UpdatePolicy { UPDATE, NOUPDATE };
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
    std::span<mgis::real> gradients;
    //! stride associate with the gradients
    const size_type gradients_stride;
    //! \brief view to the values of the thermodynamic_forces
    std::span<mgis::real> thermodynamic_forces;
    //! stride associate with the thermodynamic forces
    const size_type thermodynamic_forces_stride;
    /*!
     * \brief view to the values to the stored energy.
     */
    std::span<mgis::real> stored_energies;
    /*!
     * \brief view to the values to the dissipated energies
     */
    std::span<mgis::real> dissipated_energies;
    /*!
     * \brief material properties
     * The material properties can be uniform or not.
     * In the non uniform case, the data can be hold by the structure
     * (std::vector<real>) or simply borrow a reference
     * (std::span<mgis::real>
     * case).
     */
    std::map<std::string, FieldHolder> material_properties;
    /*! \brief mass density
     * The mass density can be uniform or not.
     * In the non uniform case, the data can be hold by the structure
     * (std::vector<real>) or simply borrow a reference
     * (std::span<mgis::real>
     * case).
     */
    std::optional<FieldHolder> mass_density;
    //! \brief view to the values of the internal state variables
    std::span<mgis::real> internal_state_variables;
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
     * (std::span<mgis::real>
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
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setMaterialProperty(
      MaterialStateManager&,
      const std::string_view&,
      const real,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
  /*!
   * \brief set the given material property
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: values
   * \param[in] s: storage mode
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setMaterialProperty(
      MaterialStateManager&,
      const std::string_view&,
      const std::span<mgis::real>&,
      const MaterialStateManager::StorageMode =
          MaterialStateManager::LOCAL_STORAGE,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
  /*!
   * \brief set the given material property
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: value
   * \param[in] p: update policy
   */
  MGIS_EXPORT [[nodiscard]] bool setMaterialProperty(
      Context&,
      MaterialStateManager&,
      const std::string_view&,
      const real,
      const MaterialStateManager::UpdatePolicy =
          MaterialStateManager::UPDATE) noexcept;
  /*!
   * \brief set the given material property
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: values
   * \param[in] s: storage mode
   * \param[in] p: update policy
   */
  MGIS_EXPORT [[nodiscard]] bool setMaterialProperty(
      Context&,
      MaterialStateManager&,
      const std::string_view&,
      const std::span<mgis::real>&,
      const MaterialStateManager::StorageMode =
          MaterialStateManager::LOCAL_STORAGE,
      const MaterialStateManager::UpdatePolicy =
          MaterialStateManager::UPDATE) noexcept;
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
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setMassDensity(
      MaterialStateManager&,
      const real,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
  /*!
   * \brief set the mass density
   * \param[out] m: material data manager
   * \param[in] v: values
   * \param[in] s: storage mode
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setMassDensity(
      MaterialStateManager&,
      const std::span<mgis::real>&,
      const MaterialStateManager::StorageMode =
          MaterialStateManager::LOCAL_STORAGE,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
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
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setExternalStateVariable(
      MaterialStateManager&,
      const std::string_view&,
      const real,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
  /*!
   * \brief set the given external state variable
   * \param[out] m: material data manager
   * \param[in] n: name
   * \param[in] v: values
   * \param[in] s: storage mode
   * \param[in] p: update policy
   */
  MGIS_EXPORT void setExternalStateVariable(
      MaterialStateManager&,
      const std::string_view&,
      const std::span<mgis::real>&,
      const MaterialStateManager::StorageMode =
          MaterialStateManager::LOCAL_STORAGE,
      const MaterialStateManager::UpdatePolicy = MaterialStateManager::UPDATE);
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
      std::span<mgis::real>,
      const mgis::behaviour::MaterialStateManager&,
      const std::string_view);

#ifdef MGIS_HAVE_HDF5

  /*!
   * \brief structure used to customize the saving of a `MaterialStateManager`
   */
  struct MaterialStateManagerSavingOptions {
    const bool allow_overwrite = true;
    const bool save_gradients = true;
    const bool save_thermodynamic_forces = true;
    const bool save_stored_energies = true;
    const bool save_dissipated_energies = true;
    const bool save_mass_densities = true;
    const bool save_material_properties = true;
    const bool save_external_state_variables = true;
  };

  /*!
   * \brief save a `MaterialStateManager` to an HDF5 group
   * \param[in] ctx: execution context
   * \param[in] g: group
   * \param[in] s: material state manager
   * \param[in] opts: options
   */
  MGIS_EXPORT [[nodiscard]] bool save(
      Context&,
      H5::Group&,
      const MaterialStateManager&,
      const MaterialStateManagerSavingOptions& = {}) noexcept;

  /*!
   * \brief structure used to customize how to restore a `MaterialStateManager`
   */
  struct MaterialStateManagerRestoreOptions {
    const bool restore_gradients = true;
    const bool restore_thermodynamic_forces = true;
    /*!
     * \brief flag stating if the stored energies shall be read
     *
     * \note this flag is ignored if the behaviour does not compute the
     * stored energy
     */
    const bool restore_stored_energies = true;
    /*!
     * \brief flag stating if the dissipated energies shall be read
     *
     * \note this flag is ignored if the behaviour does not compute the
     * dissipated energy
     */
    const bool restore_dissipated_energies = true;
    const bool restore_internal_state_variables = true;
    const bool restore_mass_densities = true;
    const bool restore_material_properties = true;
    //! \brief list of material properties that shall not be restored
    const std::vector<std::string> ignored_material_properties = {};
    const bool restore_external_state_variables = true;
    //! \brief list of external state variables that shall not be restored
    const std::vector<std::string> ignored_external_state_variables = {};
  };  // end of MaterialStateManagerRestoreOptions
  /*!
   * \brief return restore options selecting all that can be read in the given
   * group.
   *
   * \param[in] ctx: execution context
   * \param[in] g: group
   */
  MGIS_EXPORT [[nodiscard]] std::optional<MaterialStateManagerRestoreOptions>
  getGreedyMaterialStateManagerRestoreOptions(Context&,
                                              const H5::Group&) noexcept;
  /*!
   * \brief restore a `MaterialStateManager` from a HDF5 group
   *
   * \param[in] ctx: execution context
   * \param[in] g: group
   * \param[in] s: material state manager
   * \param[in] opts: options
   *
   * \note update policies are set to their values for material properties and
   * external state variables created during the restoration
   */
  MGIS_EXPORT [[nodiscard]] bool restore(
      Context&,
      MaterialStateManager&,
      const H5::Group&,
      const MaterialStateManagerRestoreOptions&) noexcept;

#endif /* MGIS_HAVE_HDF5 */

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX */
