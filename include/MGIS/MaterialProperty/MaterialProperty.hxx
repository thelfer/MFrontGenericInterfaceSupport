/*!
 * \file   include/MGIS/MaterialProperty/MaterialProperty.hxx
 * \brief
 * \author Thomas Helfer
 * \date   04/10/2022
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTY_HXX
#define LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTY_HXX

#include <string>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/MaterialProperty/MaterialPropertyFctPtr.hxx"

namespace mgis::material_property {

  /*!
   * \brief structure describing a material property
   */
  struct MGIS_EXPORT MaterialProperty {
    //! \brief constructor
    MaterialProperty();
    //! \brief move constructor
    MaterialProperty(MaterialProperty &&);
    //! \brief copy constructor
    MaterialProperty(const MaterialProperty &);
    //! \brief move assignement
    MaterialProperty &operator=(MaterialProperty &&);
    //! \brief standard assignement
    MaterialProperty &operator=(const MaterialProperty &);
    //! \brief destructor
    ~MaterialProperty();
    //! \brief name of the library providing the behaviour
    std::string library;
    //! \brief name of the material property
    std::string material_property;
    /*!
     * \brief name of the `MFront` source file used to generate the material
     * property
     */
    std::string source;
    //! \brief version of `TFEL` used to generate the material property
    std::string tfel_version;
    //! \brief unit system used by the material property
    std::string unit_system;
    //! \brief pointer to the function implementing the material property
    MaterialPropertyFctPtr fct = nullptr;
    //! \brief output of the material property
    std::string output;
    //! \brief inputs of the material property
    std::vector<std::string> inputs;
  };  // end of struct MaterialProperty

  /*!
   * \brief load the description of a material property from a library
   *
   * \param[in] l: library name
   * \param[in] b: material property name
   *
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT MaterialProperty load(const std::string &, const std::string &);

}  // end of namespace mgis::material_property

#endif /* LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTY_HXX */
