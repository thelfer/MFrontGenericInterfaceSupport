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
  constexpr BasicLinearQuadratureSpace<N>::BasicLinearQuadratureSpace(
      const size_type n) noexcept
      : nelts(n) {}

  template <size_type N>
  constexpr BasicLinearQuadratureSpace<N>::BasicLinearQuadratureSpace(
      BasicLinearQuadratureSpace&&) noexcept = default;

  template <size_type N>
  constexpr BasicLinearQuadratureSpace<N>::BasicLinearQuadratureSpace(
      const BasicLinearQuadratureSpace&) noexcept = default;

  template <size_type N>
  constexpr size_type BasicLinearQuadratureSpace<N>::size() const noexcept {
    return N * (this->nelts);
  }  // end of size

  template <size_type N>
  constexpr size_type BasicLinearQuadratureSpace<N>::getNumberOfCells()
      const noexcept {
    return this->nelts;
  }  // end of getNumberOfElements

  template <size_type N>
  constexpr size_type
  BasicLinearQuadratureSpace<N>::getNumberOfQuadraturePoints(
      const size_type) const noexcept {
    return N;
  }  // end of getNumberOfQuadraturePoints

  template <size_type N>
  constexpr typename BasicLinearQuadratureSpace<N>::DummyCellWorkspace
  BasicLinearQuadratureSpace<N>::getCellWorkspace(
      const size_type) const noexcept {
    return {};
  }  // end of getCellWorkspace

  template <size_type N>
  constexpr size_type BasicLinearQuadratureSpace<N>::getQuadraturePointOffset(
      const size_type e, const size_type i) const noexcept {
    return e * N + i;
  }  // end of getQuadraturePointOffset

  template <size_type N>
  constexpr BasicLinearQuadratureSpace<
      N>::~BasicLinearQuadratureSpace() noexcept = default;

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX */
