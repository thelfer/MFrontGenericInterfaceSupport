/*!
 * \file   MGIS/Function/BasicLinearSpace.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX
#define LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX

namespace mgis::function {

  constexpr BasicLinearSpace::BasicLinearSpace(const size_type s) noexcept
      : nelts(s) {}

  constexpr BasicLinearSpace::BasicLinearSpace(
      const BasicLinearSpace&) noexcept = default;

  constexpr BasicLinearSpace::BasicLinearSpace(BasicLinearSpace&&) noexcept =
      default;

  constexpr size_type BasicLinearSpace::size() const noexcept {
    return this->nelts;
  }

  constexpr BasicLinearSpace::~BasicLinearSpace() noexcept = default;

  constexpr size_type getSpaceSize(const BasicLinearSpace& s) noexcept {
    return s.size();
  }

  constexpr bool areEquivalent(const BasicLinearSpace& s,
                               const BasicLinearSpace& s2) noexcept {
    return s.size() == s2.size();
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX */
