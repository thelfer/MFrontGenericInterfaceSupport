/*!
 * \file   MGIS/Function/Buffer.hxx
 * \brief Definition of the `mgis::function::Buffer` type, a container used to
 * store the values of a function, and of the `makeSpan` helper functions, which
 * build a `std::span` viewing the content of such a buffer.
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

  /*!
   * \brief type used to store the values of a function
   *
   * A buffer of dynamic size is a `std::vector<real>`, whereas a buffer whose
   * size is known at compile time is a `std::array<real, Extent>`, which avoids
   * a dynamic allocation.
   *
   * \tparam Extent: number of values stored in the buffer, `dynamic_extent` if
   * the size of the buffer is only known at run time
   */
  template <size_type Extent = dynamic_extent>
  using Buffer =
      std::conditional_t<Extent == dynamic_extent,
                         std::vector<real>,
                         std::array<real, static_cast<std::size_t>(Extent)>>;

  /*!
   * \brief returns a span viewing the values stored in a fixed size buffer
   * \tparam Extent: size of the buffer
   * \param[in] b: buffer
   * \return a `std::span` viewing the values of the buffer
   * \note the extent of the returned span is only known at run time if `Extent`
   * is `std::dynamic_extent`
   */
  template <std::size_t Extent>
  auto makeSpan(const std::array<real, Extent>& b) MGIS_NOEXCEPT {
    if constexpr (Extent == std::dynamic_extent) {
      return std::span<const real>(b);
    } else {
      return std::span<const real, Extent>(b);
    }
  }  // end of makeSpan

  /*!
   * \brief returns a span viewing the values stored in a buffer of dynamic size
   * \param[in] b: buffer
   * \return a `std::span` viewing the values of the buffer
   */
  inline auto makeSpan(const std::vector<real>& b) MGIS_NOEXCEPT {
    return std::span<const real>(b);
  }  // end of makeSpan

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BUFFER_HXX */
