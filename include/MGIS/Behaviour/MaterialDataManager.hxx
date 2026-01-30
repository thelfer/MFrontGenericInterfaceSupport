/*!
 * \file   include/MGIS/Behaviour/MaterialDataManager.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_HXX
#define LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_HXX

#include <map>
#include <thread>
#include <memory>
#include <vector>
#ifdef MGIS_HAVE_HDF5
#include "MGIS/Utilities/HDF5Support.hxx"
#endif /* MGIS_HAVE_HDF5 */

#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis::behaviour {

  // forward declaration
  struct Behaviour;

  //! \brief structure in charge of handling temporary memory access.
  struct MGIS_EXPORT BehaviourIntegrationWorkSpace {
    /*!
     * \brief constructor
     * \param[in] b: behaviour
     */
    BehaviourIntegrationWorkSpace(const Behaviour&);
    //! \brief move constructor
    BehaviourIntegrationWorkSpace(BehaviourIntegrationWorkSpace&&);
    //! \brief copye constructor
    BehaviourIntegrationWorkSpace(const BehaviourIntegrationWorkSpace&);
    //! \brief move assignement
    BehaviourIntegrationWorkSpace& operator=(BehaviourIntegrationWorkSpace&&);
    //! \brief copy assignement
    BehaviourIntegrationWorkSpace& operator=(
        const BehaviourIntegrationWorkSpace&);
    //! \brief destructor
    ~BehaviourIntegrationWorkSpace();
    //! \brief a buffer to hold error messages
    std::vector<char> error_message;
    //! \brief material properties at the beginning of the time step
    std::vector<mgis::real> mps0;
    //! \brief material properties at the end of the time step
    std::vector<mgis::real> mps1;
    //! \brief external state variables at the beginning of the time step
    std::vector<mgis::real> esvs0;
    //! \brief external state variables at the end of the time step
    std::vector<mgis::real> esvs1;
    //! \brief mass density at the beginning of the time step
    mgis::real rho0;
    //! \brief mass density at the end of the time step
    mgis::real rho1;
  };  // end of struct BehaviourIntegrationWorkSpace

  /*!
   * \brief a structure in charge of holding information on how a material
   * data manager shall be initialized.
   * It may contain pointers to externally allocated data, that won't be
   * handled by the final data manager.
   * If a pointer is not initialized, the material data manager will allocate
   * and handle memory internally.
   */
  struct MaterialDataManagerInitializer {
    /*!
     * \brief view to an externally allocated memory used to store the
     * tangent operator. If empty, the material data manager will
     * initialize the required memory internally if required.
     */
    std::span<mgis::real> K;
    /*!
     * \brief view to an externally allocated memory used to store the
     * speed_of_sound. If empty, the material data manager will
     * initialize the required memory internally if required.
     */
    std::span<mgis::real> speed_of_sound;
    /*!
     * \brief object used to initalize the state manager associated with the
     * beginning of the time step.
     */
    MaterialStateManagerInitializer s0;
    /*!
     * \brief object used to initalize the state manager associated with the
     * end of the time step.
     */
    MaterialStateManagerInitializer s1;
  };  // end of MaterialDataManagerInitializer

  /*!
   * \brief structure in charge of handling the data associated with a
   * material in an optimized way. Here, the "material" is defined by a
   * behaviour and a number of integration points.
   *
   * The following design choices were made:
   * - The material properties and the external state variables are treated
   *   individually. They can be uniform or spatially variable.
   * - The internal state variables are treated as a block.
   */
  struct MGIS_EXPORT MaterialDataManager {
    /*!
     * \brief main constructor
     * \param[in] behaviour: behaviour
     * \param[in] s: number of integration points
     */
    MaterialDataManager(const Behaviour&, const size_type);
    /*!
     * \brief main constructor
     * \param[in] behaviour: behaviour
     * \param[in] s: number of integration points
     * \param[in] i: initializer
     */
    MaterialDataManager(const Behaviour&,
                        const size_type,
                        const MaterialDataManagerInitializer&);
    /*!
     * \brief set if the `MaterialDataManager` must take care of thread-safety.
     * This flag is mostly used in members functions allocating memory.
     * \param[in] bv: boolean
     */
    void setThreadSafe(const bool);
    /*!
     * \brief allocate the memory associated with the tangent operator blocks if
     * required.
     *
     * This method is useless if the memory associated with the tangent operator
     * blocks had previously been allocated or assigned to external memory (see
     * the `MaterialDataManagerInitializer` structure).
     *
     * \note This method is thread-safe if the `thread_safe` is `true`.
     * In this case, the memory allocation is guarded by a mutex.
     * See the `setThreadSafe` method for details
     */
    void allocateArrayOfTangentOperatorBlocks();
    /*!
     * \brief use an externally allocated memory to store the tangent operator
     * blocks.
     *
     * \param[in] m: memory view
     *
     * \note this method calls `releaseArrayOfTangentOperatorBlocks` before
     * allocating the memory.
     */
    void useExternalArrayOfTangentOperatorBlocks(std::span<real>);
    /*!
     * \brief release the memory associated with the tangent operator blocks.
     *
     * If the memory associated with the tangent operator blocks was handled
     * internally, this memory is freed.
     * If an external memory buffer was used, reference to this buffer is
     * removed.
     */
    void releaseArrayOfTangentOperatorBlocks();
    /*!
     * \brief allocate the memory associated with the speed of sound if
     * required.
     *
     * This method is useless if the memory associated with the speed of sound
     * had previously been allocated or assigned to external memory (see the
     * `MaterialDataManagerInitializer` structure).
     *
     * \note This method is thread-safe if the `thread_safe` is `true`.
     * In this case, the memory allocation is guarded by a mutex.
     * See the `setThreadSafe` method for details
     */
    void allocateArrayOfSpeedOfSounds();
    /*!
     * \brief use an externally allocated memory to store the tangent operator
     * blocks.
     *
     * \param[in] m: memory view
     *
     * \note this method calls `releaseArrayOfSpeedOfSounds` before
     * allocating the memory.
     */
    void useExternalArrayOfSpeedOfSounds(std::span<real>);
    /*!
     * \brief release the memory associated with the tangent operator blocks.
     *
     * If the memory associated with the tangent operator blocks was handled
     * internally, this memory is freed.
     * If an external memory buffer was used, reference to this buffer is
     * removed.
     */
    void releaseArrayOfSpeedOfSounds();
    /*!
     * \brief return a workspace associated with the given behaviour.
     *
     * \note This method returns a object per thread if the `thread_safe` member
     * is `true`.
     */
    BehaviourIntegrationWorkSpace& getBehaviourIntegrationWorkSpace();
    /*!
     * \brief clear behaviour integration workspaces.
     *
     * \note This is only interesting if you keep allocating different thread
     * pools.
     */
    void releaseBehaviourIntegrationWorkspaces();
    //! \brief destructor
    ~MaterialDataManager();
    //! \brief state at the beginning of the time step
    MaterialStateManager s0;
    //! \brief state at the end of the time step
    MaterialStateManager s1;
    //! \brief view of the stiffness matrices, if required.
    std::span<real> K;
    /*!
     * \brief proposed time step increment increase factor
     *
     * The calling solver shall set a suitable value on input
     * depending on its policy **before** each call to integrate.
     *
     * For instance, if the solver want to limit the increase to 20% at most, it
     * shall set it to 1.2. But setting it to 1, the solver won't allow the
     * behaviour to request an increase of the time step.
     */
    real rdt = 1;
    //! \brief view on the speed of sound.
    std::span<real> speed_of_sound;
    //! \brief number of integration points
    const size_type n;
    /*!
     * \brief the size of the stiffness matrix for one integration point (the
     * size of K is K_stride times the number of integration points)
     */
    const size_type K_stride;
    //! \brief underlying behaviour
    const Behaviour& b;

   private:
    //! move constructor
    MaterialDataManager(MaterialDataManager&&) = delete;
    //! copy constructor
    MaterialDataManager(const MaterialDataManager&) = delete;
    //! move assignement
    MaterialDataManager& operator=(MaterialDataManager&&) = delete;
    //! copy assignement
    MaterialDataManager& operator=(const MaterialDataManager&) = delete;
    //! \brief values of the stiffness matrices, if hold internally.
    std::vector<real> K_values;
    //! \brief values of the speed of sound, if hold internally.
    std::vector<real> speed_of_sound_values;
    //! \brief integration workspace for individual threads.
    std::map<std::thread::id, std::unique_ptr<BehaviourIntegrationWorkSpace>>
        iwks;
    //! \brief a pointer to an integration workspace
    std::unique_ptr<BehaviourIntegrationWorkSpace> iwk;
    //! \brief boolean stating if thread safety must be unsured
    bool thread_safe = true;
  };  // end of struct MaterialDataManager

  /*!
   * \brief update the behaviour data by:
   * - setting s1 equal to s0
   * - filling the stiffness matrix with 0
   * \param[in,out] m: material data manager
   */
  MGIS_EXPORT void update(MaterialDataManager&);
  /*!
   * \brief revert the behaviour data by:
   * - setting s1 equal to s0
   * - filling the stiffness matrix with 0
   * \param[in,out] m: material data manager
   */
  MGIS_EXPORT void revert(MaterialDataManager&);

  /*!
   * \return an array containing the results of a post-processing.
   * \param[in] m: material data manager
   * \param[in] n: name of the post-processing
   */
  MGIS_EXPORT std::vector<mgis::real> allocatePostProcessingVariables(
      const MaterialDataManager&, const std::string_view);

#ifdef MGIS_HAVE_HDF5

  /*!
   * \brief structure used to customize the saving of a `MaterialDataManager`
   */
  struct MaterialDataManagerSavingOptions
      : public MaterialStateManagerSavingOptions {};

  /*!
   * \brief save a `MaterialDataManager` to a file
   * \param[in] ctx: execution context
   * \param[in] g: group
   * \param[in] m: material data manager
   * \param[in] opts: options
   */
  MGIS_EXPORT [[nodiscard]] bool save(
      Context&,
      H5::Group&,
      const MaterialDataManager&,
      const MaterialDataManagerSavingOptions& = {}) noexcept;

#endif /* MGIS_HAVE_HDF5 */

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_HXX */
