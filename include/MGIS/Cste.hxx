/*!
 * \file   Cste.hxx
 * \brief
 * \author Thomas Helfer
 * \date   03/09/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_CSTE_HXX
#define LIB_MGIS_CSTE_HXX

#include <limits>
#include "MGIS/Config.hxx"

namespace mgis {

  namespace internals {

    /*!
     * \brief compute the absolute value of a number using a constexpr
     * function.
     */
    template <typename T>
    constexpr T abs(const T v) {
      return (v < 0) ? -v : v;
    }

    template <typename T>
    constexpr T sqrt_helper(const T x, const T g) {
      return abs(g - x / g) < 2 * std::numeric_limits<T>::epsilon()
                 ? g
                 : sqrt_helper(x, (g + x / g) / T{2});
    }

    /*!
     * \brief compute the square root of a number using the Heron
     * algorithm using a constexpr function.
     *
     * This function is meant to compute some usefull constant such as
     * sqrt(2) and sqrt(3) at compile-time.
     * \tparam T : numerical type
     * \pre T must be a floatting point type
     * \pre std::numeric_limits<T>::epsilon must be defined
     */
    template <typename T>
    constexpr T sqrt(const T v) {
      return sqrt_helper(v, T{1});
    }

  }  // end of namespace internals

  struct Cste {
#ifndef _MSC_VER
    static constexpr real sqrt2 = internals::sqrt(real{2});
    static constexpr real isqrt2 = 1 / internals::sqrt(real{2});
    static constexpr real sqrt3 = internals::sqrt(real{3});
    static constexpr real isqrt3 = 1 / internals::sqrt(real{3});
#else
    static constexpr const real sqrt2 = 1.41421356237309504880;
    static constexpr const real isqrt2 = 0.70710678118654752440;
    static constexpr const real sqrt3 =
        1.7320508075688772935274463415058723669428052538103806280;
    static constexpr const real isqrt3 = 0.57735026919;
#endif
  };

}  // end of namespace mgis

#endif /* LIB_MGIS_CSTE_HXX */
