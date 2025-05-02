/*!
 * \file   MGIS/Buffer.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   30/04/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_BUFFER_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_BUFFER_HXX

#include <span>
#include <array>
#include <vector>
#include <type_traits>
#include <MGIS/Config.hxx>

namespace mgis::quadrature_function {

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

}  // namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_BUFFER_HXX */
