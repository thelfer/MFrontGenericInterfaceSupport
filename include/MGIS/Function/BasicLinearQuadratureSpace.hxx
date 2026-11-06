/*!
 * \file   MGIS/Function/BasicLinearQuadratureSpace.hxx
 * \brief Definition of `mgis::function::BasicLinearQuadratureSpace`, a minimal
 * implementation of a quadrature space whose cells all have the same number of
 * integration points, and of the functions associated with this space.
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

#ifndef LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX
#define LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a minimal implementation of a linear quadrature space
   *
   * The elements of such a space are the quadrature points of all the cells of
   * the space. Given the index `e` of a cell and the index `i` of one of its
   * quadrature points, the offset of the corresponding element is `e * N + i`,
   * which is why no workspace is required to access the cells of the space.
   *
   * \tparam N: number of integration points per cell, which must be greater
   * than zero
   */
  template <size_type N>
  requires(N > 0) struct BasicLinearQuadratureSpace {
    /*!
     * \brief an empty cell workspace, as the cells of this space do not require
     * any information to be accessed
     */
    struct DummyCellWorkspace {};
    /*!
     * \brief constructor
     * \param[in] n: number of cells of the space
     */
    constexpr BasicLinearQuadratureSpace(const size_type) MGIS_NOEXCEPT;
    //! \brief move constructor
    constexpr BasicLinearQuadratureSpace(BasicLinearQuadratureSpace&&) noexcept;
    //! \brief copy constructor
    constexpr BasicLinearQuadratureSpace(const BasicLinearQuadratureSpace&)
        MGIS_NOEXCEPT;
    /*!
     * \return the number of quadrature points of the space, i.e. the number of
     * cells multiplied by `N`
     */
    [[nodiscard]] constexpr size_type size() const MGIS_NOEXCEPT;
    /*!
     * \return the number of cells in the quadrature space
     */
    [[nodiscard]] constexpr size_type getNumberOfCells() const MGIS_NOEXCEPT;
    /*!
     * \brief return the number of quadrature points of a given cell
     * \param[in] e: cell index
     * \return the number of quadrature points of the given cell, which does
     * not depend on the cell index `e` and is always `N`
     */
    [[nodiscard]] constexpr size_type getNumberOfQuadraturePoints(
        const size_type) const MGIS_NOEXCEPT;
    /*!
     * \brief return the workspace of a given cell
     * \param[in] e: cell index
     * \return the workspace of the cell of index `e`, which is always an
     * empty `DummyCellWorkspace`
     */
    [[nodiscard]] constexpr DummyCellWorkspace getCellWorkspace(
        const size_type) const MGIS_NOEXCEPT;
    /*!
     * \brief return the offset of a given quadrature point in the space
     * \param[in] e: cell index
     * \param[in] i: quadrature point index of the cell of index `e`
     * \return the offset in the space of the quadrature point of index `i` of
     * the cell of index `e`, that is `e * N + i`
     */
    [[nodiscard]] constexpr size_type getQuadraturePointOffset(
        const size_type, const size_type) const MGIS_NOEXCEPT;
    //! \brief destructor
    constexpr ~BasicLinearQuadratureSpace() noexcept;

   private:
    //! \brief number of cells of the space
    size_type ncells;
  };

  /*!
   * \brief traits defining the properties of a `BasicLinearQuadratureSpace`
   *
   * \tparam N: number of integration points per cell
   */
  template <size_type N>
  struct SpaceTraits<BasicLinearQuadratureSpace<N>> {
    /*!
     * \brief the elements of the space are addressed by an index varying from 0
     * to the size of the space minus one
     */
    static constexpr bool linear_element_indexing = true;
    /*!
     * \brief the cells of the space are addressed by an index varying from 0 to
     * the number of cells minus one
     */
    static constexpr bool linear_cell_indexing = true;
    //! \brief type used to store the size of the space
    using size_type = mgis::size_type;
    //! \brief type of the index of an element of the space
    using element_index_type = mgis::size_type;
    //! \brief type of the index of a cell of the space
    using cell_index_type = mgis::size_type;
    //! \brief type of the index of a quadrature point of a cell
    using quadrature_point_index_type = mgis::size_type;
  };

  /*!
   * \brief returns the size of the space, i.e. the number of quadrature points
   * it contains
   * \tparam N: number of integration points per cell
   * \param[in] s: space
   * \return the number of elements of the space
   */
  template <size_type N>
  constexpr size_type getSpaceSize(const BasicLinearQuadratureSpace<N>&)
      MGIS_NOEXCEPT;

  /*!
   * \brief returns the number of cells of the space
   * \tparam N: number of integration points per cell
   * \param[in] s: space
   * \return the number of cells of the space
   */
  template <size_type N>
  constexpr size_type getNumberOfCells(const BasicLinearQuadratureSpace<N>&)
      MGIS_NOEXCEPT;

  /*!
   * \brief returns the number of quadrature points of a given cell, which is
   * always `N`
   * \tparam N: number of integration points per cell
   * \param[in] s: space
   * \param[in] e: cell index
   * \return the number of quadrature points of the cell of index `e`
   */
  template <size_type N>
  constexpr size_type getNumberOfQuadraturePoints(
      const BasicLinearQuadratureSpace<N>&, const size_type) MGIS_NOEXCEPT;

  /*!
   * \brief returns the offset in the space of a given quadrature point
   * \tparam N: number of integration points per cell
   * \param[in] s: space
   * \param[in] e: cell index
   * \param[in] i: quadrature point index of the cell of index `e`
   * \return the offset `e * N + i` of the quadrature point
   */
  template <size_type N>
  constexpr size_type getQuadraturePointOffset(
      const BasicLinearQuadratureSpace<N>&,
      const size_type,
      const size_type) MGIS_NOEXCEPT;

  /*!
   * \brief checks if two spaces are equivalent, i.e. if they have the same
   * number of quadrature points
   * \tparam N: number of integration points per cell
   * \param[in] s: first space
   * \param[in] s2: second space
   * \return true if the two spaces have the same size
   */
  template <size_type N>
  [[nodiscard]] constexpr bool areEquivalent(
      const BasicLinearQuadratureSpace<N>&,
      const BasicLinearQuadratureSpace<N>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX */
