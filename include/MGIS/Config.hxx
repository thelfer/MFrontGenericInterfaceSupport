/*!
 * \file   include/MGIS/Config.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_CONFIG_HXX
#define LIB_MGIS_CONFIG_HXX

#include <limits>
#include "MGIS/Config-c.h"

namespace mgis {

  namespace attributes {

    /*!
     * \brief an attribute use to indicate that a method or a function as being
     * unsafe without precautions
     */
    struct UnsafeAttribute {};
    /*!
     * \brief an attribute use to indicate that a method may throw or not
     */
    template <bool>
    struct ThrowingAttribute {};

  }  // namespace attributes

  //
  inline constexpr auto unsafe = attributes::UnsafeAttribute{};
  //
  inline constexpr auto throwing = attributes::ThrowingAttribute<true>{};
  //
  inline constexpr auto not_throwing = attributes::ThrowingAttribute<false>{};

  //! \brief a simple alias to the the default indexing type used by mgis
  using size_type = mgis_size_type;

  //! \brief alias to the numeric type used
  using real = mgis_real;

  //! \brief a constant whose role is similar to std::dynamic_extent
  inline constexpr size_type dynamic_extent =
      std::numeric_limits<size_type>::max();

}  // end of namespace mgis

#endif /* LIB_MGIS_CONFIG_HXX */
