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
      //! \brief values of the thermodynamic_forces
      std::vector<real> thermodynamic_forces;
      //! \brief material properties
      std::map<std::string, mgis::variant<real, std::vector<real>>>
          material_properties;
      //! \brief values of the internal state variables
      std::vector<real> internal_state_variables;
      //! \brief values of the external state variables
      std::map<std::string, mgis::variant<real, std::vector<real>>>
          external_state_variables;
      //! underlying behaviour
      const Behaviour& b;
    };  // end of struct MaterialStateManager

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALSTATEMANAGER_HXX */
