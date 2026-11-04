/*!
 * \file   MGIS/Function/BasicLinearSpace.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX
#define LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX

namespace mgis::function {

  constexpr BasicLinearSpace::BasicLinearSpace(const size_type s) MGIS_NOEXCEPT
      : nelts(s) {}

  constexpr BasicLinearSpace::BasicLinearSpace(
      const BasicLinearSpace&) MGIS_NOEXCEPT = default;

  constexpr BasicLinearSpace::BasicLinearSpace(BasicLinearSpace&&) noexcept =
      default;

  constexpr size_type BasicLinearSpace::size() const MGIS_NOEXCEPT {
    return this->nelts;
  }

  constexpr BasicLinearSpace::~BasicLinearSpace() noexcept = default;

  constexpr size_type getSpaceSize(const BasicLinearSpace& s) MGIS_NOEXCEPT {
    return s.size();
  }

  constexpr bool areEquivalent(const BasicLinearSpace& s,
                               const BasicLinearSpace& s2) MGIS_NOEXCEPT {
    return s.size() == s2.size();
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BASICLINEARSPACE_IXX */
