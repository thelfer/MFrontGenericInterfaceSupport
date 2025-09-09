/*!
 * \file   MGIS/Function/BasicLinearQuadratureSpace.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX
#define LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX

namespace mgis::function {

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::BasicLinearQuadratureSpace(const size_type n) noexcept
      : ncells(n) {}

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::BasicLinearQuadratureSpace(BasicLinearQuadratureSpace&&) noexcept =
      default;

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<N>::
      BasicLinearQuadratureSpace(const BasicLinearQuadratureSpace&) noexcept =
          default;

  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::size() const noexcept {
    return N * (this->ncells);
  }  // end of size

  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getNumberOfCells() const noexcept {
    return this->ncells;
  }  // end of getNumberOfElements

  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getNumberOfQuadraturePoints(
          const size_type) const noexcept {
    return N;
  }  // end of getNumberOfQuadraturePoints

  template <size_type N>
  requires(N > 0) constexpr
      typename BasicLinearQuadratureSpace<N>::DummyCellWorkspace
      BasicLinearQuadratureSpace<N>::getCellWorkspace(
          const size_type) const noexcept {
    return {};
  }  // end of getCellWorkspace

  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getQuadraturePointOffset(
          const size_type e, const size_type i) const noexcept {
    return e * N + i;
  }  // end of getQuadraturePointOffset

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::~BasicLinearQuadratureSpace() noexcept = default;

  template <size_type N>
  constexpr size_type getSpaceSize(
      const BasicLinearQuadratureSpace<N>& s) noexcept {
    return s.size();
  }

  template <size_type N>
  constexpr size_type getNumberOfCells(
      const BasicLinearQuadratureSpace<N>& s) noexcept {
    return s.getNumberOfCells();
  }  // end of getNumberOfCells

  template <size_type N>
  constexpr size_type getNumberOfQuadraturePoints(
      const BasicLinearQuadratureSpace<N>& s, const size_type e) noexcept {
    return s.getNumberOfQuadraturePoints(e);
  }  // end of getNumberOfQuadraturePoints

  template <size_type N>
  constexpr size_type getQuadraturePointOffset(
      const BasicLinearQuadratureSpace<N>& s,
      const size_type e,
      const size_type i) noexcept {
    return s.getQuadraturePointOffset(e, i);
  }

  template <size_type N>
  constexpr bool areEquivalent(
      const BasicLinearQuadratureSpace<N>& s,
      const BasicLinearQuadratureSpace<N>& s2) noexcept {
    return s.size() == s2.size();
  }

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX */
