/*!
 * \file   MGIS/Function/BasicLinearQuadratureSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX
#define LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a minimal implementation of linear quadrature space
   * \param N: number of integration points per cells
   */
  template <size_type N>
  requires(N > 0) struct MGIS_EXPORT BasicLinearQuadratureSpace {
    //! \brief an empty element workspace
    struct DummyCellWorkspace {};
    /*!
     * \brief constructor
     * \param[in] n: number of elements
     */
    constexpr BasicLinearQuadratureSpace(const size_type) noexcept;
    constexpr BasicLinearQuadratureSpace(BasicLinearQuadratureSpace&&) noexcept;
    constexpr BasicLinearQuadratureSpace(
        const BasicLinearQuadratureSpace&) noexcept;
    //! \return the number of quadrature points
    constexpr size_type size() const noexcept;
    //! \return the number of cells in the quadrature space
    constexpr size_type getNumberOfCells() const noexcept;
    /*!
     * \return the number quadrature points for the given element
     * \param[in] e: cell index
     */
    constexpr size_type getNumberOfQuadraturePoints(
        const size_type) const noexcept;
    /*!
     * \return the workspace of the element
     * \param[in] e: cell index
     */
    constexpr DummyCellWorkspace getCellWorkspace(
        const size_type) const noexcept;
    /*!
     * \return the offset associated with a quadrature space
     * \param[in] wk: element workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     */
    constexpr size_type getQuadraturePointOffset(
        const size_type, const size_type) const noexcept;
    //! \brief destructor
    constexpr ~BasicLinearQuadratureSpace() noexcept;

   private:
    //! \brief number of cells of the the space
    size_type ncells;
  };

  template <size_type N>
  struct SpaceTraits<BasicLinearQuadratureSpace<N>> {
    static constexpr bool linear_element_indexing = true;
    static constexpr bool linear_cell_indexing = true;
    using size_type = mgis::size_type;
    using element_index_type = mgis::size_type;
    using cell_index_type = mgis::size_type;
    using quadrature_point_index_type = mgis::size_type;
  };

  template <size_type N>
  constexpr size_type getSpaceSize(
      const BasicLinearQuadratureSpace<N>&) noexcept;

  template <size_type N>
  constexpr size_type getNumberOfCells(
      const BasicLinearQuadratureSpace<N>&) noexcept;

  template <size_type N>
  constexpr size_type getNumberOfQuadraturePoints(
      const BasicLinearQuadratureSpace<N>&, const size_type) noexcept;

  template <size_type N>
  constexpr size_type getQuadraturePointOffset(
      const BasicLinearQuadratureSpace<N>&,
      const size_type,
      const size_type) noexcept;

  template <size_type N>
  [[nodiscard]] constexpr bool areEquivalent(
      const BasicLinearQuadratureSpace<N>&,
      const BasicLinearQuadratureSpace<N>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX */
