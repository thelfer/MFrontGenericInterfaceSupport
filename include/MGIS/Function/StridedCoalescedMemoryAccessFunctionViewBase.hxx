/*!
 * \file   MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.hxx
 * \brief
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
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) struct StridedCoalescedMemoryAccessFunctionViewBase
      : private PreconditionsChecker<
            StridedCoalescedMemoryAccessFunctionViewBase<Space,
                                                         N,
                                                         is_mutable>> {
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] space: space
     * \param[in] values: values
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space& space,
        std::span<const real>) requires(!is_mutable);
    /*!
     * \param[in] space: space
     * \param[in] values: values
     */
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const Space&, std::span<const real>) requires(!is_mutable);
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] space: space
     * \param[in] values: values
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&, const Space&, std::span<real>);
    /*!
     * \param[in] space: space
     * \param[in] values: values
     */
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(const Space&,
                                                           std::span<real>);
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr std::array<real, N> getValues(const size_type) const
        requires(LinearElementSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr std::array<real, N> getValues(const size_type,
                                                          const size_type) const
        requires(LinearQuadratureSpaceConcept<Space>);
    //! \return the underlying quadrature space
    [[nodiscard]] constexpr const Space& getSpace() const noexcept;

   protected:
    /*!
     * \param[in] space: space
     * \param[in] values: values
     */
    template <bool doPreconditionsCheck>
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<const real>) requires(!is_mutable);
    /*!
     * \param[in] space: space
     * \param[in] values: values
     */
    template <bool doPreconditionsCheck>
    constexpr StridedCoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<real>);
    //! \brief underlying discretization space
    const Space space;
    //! \brief pointer to the first component
    const std::conditional_t<is_mutable, real*, const real*> data_pointer;
  };

}  // namespace mgis::function

#include "MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.ixx"

#endif /* LIB_MGIS_FUNCTION_STRIDEDCOALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX \
        */
