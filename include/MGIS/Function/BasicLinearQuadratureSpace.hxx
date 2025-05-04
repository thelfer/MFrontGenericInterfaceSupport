/*!
 * \file   MGIS/Function/BasicLinearQuadratureSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX
#define LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function/Space.hxx"

namespace mgis::function {

  /*!
   * \brief a minimal implementation of linear quadrature space
   * \param N: number of integration points per element
   */
  template <size_type N>
  requires(N > 0) struct MGIS_EXPORT BasicLinearQuadratureSpace {
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
    //! \return the number of elements
    constexpr size_type getNumberOfElements() const noexcept;
    /*!
     * \return the number quadrature points for the given element
     * \param[in] e: element index
     */
    constexpr size_type getNumberOfQuadraturePoints(
        const size_type) const noexcept;
    //! \brief destructor
    constexpr ~BasicLinearQuadratureSpace() noexcept;

   private:
    //! \brief number of elements of the the space
    size_type nelts;
  };

  template <size_type N>
  struct SpaceTraits<BasicLinearQuadratureSpace<N>> {
    using size_type = mgis::size_type;
    using element_index_type = mgis::size_type;
    using quadrature_point_index_type = mgis::size_type;
  };

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearQuadratureSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARQUADRATURESPACE_HXX */
