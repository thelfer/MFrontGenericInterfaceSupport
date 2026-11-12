/*!
 * \file   MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.hxx
 * \brief File which defines a base class for views whose components are
 * stored contiguously in a single memory region, in a non interleaved manner.
 * \author Thomas Helfer
 * \date   17/01/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_STRIDEDCOALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX
#define LIB_MGIS_FUNCTION_STRIDEDCOALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX

#include <span>
#include <array>
#include "MGIS/Config.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  /*!
   * \brief a class meant to describe a function, each component of which is
   * stored contiguously in non interleaved manner using the following scheme:
   *
   * ~~~~
   * | <------- Component 1 ---------> |....| <----- Component Nc ---------> |
   * +-------++-------++------++-------+----+-------++------++------++-------+
   * | Elt 1 || Elt 2 || .... || Elt N |....| Elt 1 ||Elt 2 || .... || Elt N |
   * +-------++-------++------++-------+----+-------++------++------++-------+
   * ~~~~
   *
   * Contrary to the `CoalescedMemoryAccessFunctionViewBase` class, the values
   * of all the components are stored in a single contiguous memory region:
   * the values of the `k`-th component start at the offset
   * `k * getSpaceSize(space)`, i.e. the stride between the values of two
   * successive components at a given integration point is equal to the size
   * of the space. This allows the class to store only the space and a pointer
   * to the first value.
   *
   * \tparam Space: type of the space on which the view is defined
   * \tparam N: number of components, which must be strictly positive
   * \tparam is_mutable: boolean stating if the data pointer is mutable
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) struct StridedCoalescedMemoryAccessFunctionViewBase
      : private PreconditionsChecker<
            StridedCoalescedMemoryAccessFunctionViewBase<Space,
                                                         N,
                                                         is_mutable>> {
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] space: space on which the view is defined
     * \param[in] values: constant values
     * \return true if the number of values is equal to
     * `N * getSpaceSize(space)`
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space& space,
        std::span<const real>) requires(!is_mutable);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] space: space on which the view is defined
     * \param[in] values: constant values
     */
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const Space&, std::span<const real>) requires(!is_mutable);
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] space: space on which the view is defined
     * \param[in] values: values
     * \return true if the number of values is equal to
     * `N * getSpaceSize(space)`
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&, const Space&, std::span<real>);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] space: space on which the view is defined
     * \param[in] values: values
     */
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(const Space&,
                                                           std::span<real>);
    /*!
     * \return the values of all the components at the given integration point
     * \param[in] i: offset associated with the integration point
     * \note the values are copied in a `std::array` as the values of the
     * different components at a given integration point are not contiguous in
     * memory
     */
    [[nodiscard]] constexpr std::array<real, N> getValues(const size_type) const
        requires(LinearElementSpaceConcept<Space>);
    /*!
     * \return the values of all the components at the given integration point
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note the values are copied in a `std::array` as the values of the
     * different components at a given integration point are not contiguous in
     * memory
     */
    [[nodiscard]] constexpr std::array<real, N> getValues(const size_type,
                                                          const size_type) const
        requires(LinearQuadratureSpaceConcept<Space>);
    /*!
     * \brief return the underlying space
     * \return the space on which the view is defined
     */
    [[nodiscard]] constexpr const Space& getSpace() const MGIS_NOEXCEPT;

   protected:
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] space: space on which the view is defined
     * \param[in] values: constant values
     */
    template <bool doPreconditionsCheck>
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<const real>) requires(!is_mutable);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] space: space on which the view is defined
     * \param[in] values: values
     */
    template <bool doPreconditionsCheck>
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<real>);
    //! \brief space on which the view is defined
    const Space space;
    //! \brief pointer to the first value, i.e. the value of the first
    //! component of the first point of the space
    const std::conditional_t<is_mutable, real*, const real*> data_pointer;
  };

}  // namespace mgis::function

#include "MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.ixx"

#endif /* LIB_MGIS_FUNCTION_STRIDEDCOALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX \
        */
