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

#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Variant.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis {

  namespace behaviour {

    // forward declaration
    struct Behaviour;

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
       * initialize the required memory internally.
       */
      mgis::span<mgis::real> K;
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
      //! destructor
      ~MaterialDataManager();
      //! \brief state at the beginning of the time step
      MaterialStateManager s0;
      //! \brief state at the end of the time step
      MaterialStateManager s1;
      //! \brief view of the stiffness matrices.
      mgis::span<real> K;
      //! \brief number of integration points
      const size_type n;
      /*!
       * \brief the size of the stiffness matrix for one integration point (the
       * size of K is K_stride times the number of integration points)
       */
      const size_type K_stride;
      //! underlying behaviour
      const Behaviour& b;

     private:
      //! \brief values of the stiffness matrices, if hold internally.
      std::vector<real> K_values;
      //! move constructor
      MaterialDataManager(MaterialDataManager&&) = delete;
      //! copy constructor
      MaterialDataManager(const MaterialDataManager&) = delete;
      //! move assignement
      MaterialDataManager& operator=(MaterialDataManager&&) = delete;
      //! copy assignement
      MaterialDataManager& operator=(const MaterialDataManager&) = delete;
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

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALDATAMANAGER_HXX */
