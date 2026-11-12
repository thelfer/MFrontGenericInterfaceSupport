/*!
 * \file   MGIS/Function/BasicLinearQuadratureSpace.ixx
 * \brief Implementation of the `mgis::function::BasicLinearQuadratureSpace`
 * structure and of the functions associated with this space.
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX
#define LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX

namespace mgis::function {

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::BasicLinearQuadratureSpace(const size_type n) MGIS_NOEXCEPT
      : ncells(n) {}

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::BasicLinearQuadratureSpace(BasicLinearQuadratureSpace&&) noexcept =
      default;

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<N>::
      BasicLinearQuadratureSpace(const BasicLinearQuadratureSpace&) MGIS_NOEXCEPT =
          default;

  /*!
   * \return the total number of quadrature points of the space, which is the
   * number of cells multiplied by the number of quadrature points per cell
   */
  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::size() const MGIS_NOEXCEPT {
    return N * (this->ncells);
  }  // end of size

  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getNumberOfCells() const MGIS_NOEXCEPT {
    return this->ncells;
  }  // end of getNumberOfElements

  /*!
   * \return the number of quadrature points of the considered cell, which does
   * not depend on the cell since all the cells have `N` quadrature points
   */
  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getNumberOfQuadraturePoints(
          const size_type) const MGIS_NOEXCEPT {
    return N;
  }  // end of getNumberOfQuadraturePoints

  /*!
   * \return an empty workspace: the quadrature points of a cell of this space
   * are located from the index of this cell only
   */
  template <size_type N>
  requires(N > 0) constexpr
      typename BasicLinearQuadratureSpace<N>::DummyCellWorkspace
      BasicLinearQuadratureSpace<N>::getCellWorkspace(
          const size_type) const MGIS_NOEXCEPT {
    return {};
  }  // end of getCellWorkspace

  /*!
   * \return the offset of the quadrature point: the cells are stored one after
   * the other and all have `N` quadrature points
   */
  template <size_type N>
  requires(N > 0) constexpr size_type
      BasicLinearQuadratureSpace<N>::getQuadraturePointOffset(
          const size_type e, const size_type i) const MGIS_NOEXCEPT {
    return e * N + i;
  }  // end of getQuadraturePointOffset

  template <size_type N>
  requires(N > 0) constexpr BasicLinearQuadratureSpace<
      N>::~BasicLinearQuadratureSpace() noexcept = default;

  //! \brief forwards the call to the `size` member function
  template <size_type N>
  constexpr size_type getSpaceSize(
      const BasicLinearQuadratureSpace<N>& s) MGIS_NOEXCEPT {
    return s.size();
  }

  //! \brief forwards the call to the `getNumberOfCells` member function
  template <size_type N>
  constexpr size_type getNumberOfCells(
      const BasicLinearQuadratureSpace<N>& s) MGIS_NOEXCEPT {
    return s.getNumberOfCells();
  }  // end of getNumberOfCells

  //! \brief forwards the call to the `getNumberOfQuadraturePoints` member
  /// function
  template <size_type N>
  constexpr size_type getNumberOfQuadraturePoints(
      const BasicLinearQuadratureSpace<N>& s, const size_type e) MGIS_NOEXCEPT {
    return s.getNumberOfQuadraturePoints(e);
  }  // end of getNumberOfQuadraturePoints

  //! \brief forwards the call to the `getQuadraturePointOffset` member function
  template <size_type N>
  constexpr size_type getQuadraturePointOffset(
      const BasicLinearQuadratureSpace<N>& s,
      const size_type e,
      const size_type i) MGIS_NOEXCEPT {
    return s.getQuadraturePointOffset(e, i);
  }

  // two spaces are declared equivalent when they have the same number of
  // quadrature points, no identity of the spaces is required
  template <size_type N>
  constexpr bool areEquivalent(
      const BasicLinearQuadratureSpace<N>& s,
      const BasicLinearQuadratureSpace<N>& s2) MGIS_NOEXCEPT {
    return s.size() == s2.size();
  }

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_IXX */
