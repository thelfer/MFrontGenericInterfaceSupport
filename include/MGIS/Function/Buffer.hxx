/*!
 * \file   MGIS/Buffer.hxx
 * \brief
 * \author Thomas Helfer
 * \date   30/04/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_BUFFER_HXX
#define LIB_MGIS_FUNCTION_BUFFER_HXX

#include <span>
#include <array>
#include <vector>
#include <type_traits>
#include <MGIS/Config.hxx>

namespace mgis::function {

  template <size_type Extent = dynamic_extent>
  using Buffer =
      std::conditional_t<Extent == dynamic_extent,
                         std::vector<real>,
                         std::array<real, static_cast<std::size_t>(Extent)>>;

  template <std::size_t Extent>
  auto makeSpan(const std::array<real, Extent>& b) noexcept {
    if constexpr (Extent == std::dynamic_extent) {
      return std::span<const real>(b);
    } else {
      return std::span<const real, Extent>(b);
    }
  }  // end of makeSpan

  inline auto makeSpan(const std::vector<real>& b) noexcept {
    return std::span<const real>(b);
  }  // end of makeSpan

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BUFFER_HXX */
