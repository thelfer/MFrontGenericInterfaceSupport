/*!
 * \file   MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/10/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX
#define LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX

#include <span>
#include <array>
#include "MGIS/Config.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable = true>
  requires(N > 0) struct CoalescedMemoryAccessFunctionViewBase
      : private PreconditionsChecker<
            CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>> {
    //
    using MutableValues = std::array<real*, N>;
    //
    using ConstValues = std::array<const real*, N>;

    //! \brief type of the function view associated with a single component
    using ScalarComponentFunctionView =
        FunctionView<Space,
                     FunctionDataLayoutDescription{.data_size = 1,
                                                   .data_stride = 1},
                     is_mutable>;
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] components: components
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \param[in] components: components
     */
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \param[in] components: components
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr MutableValues
    getValuesPointers(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr MutableValues
    getValuesPointers(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr ConstValues getValuesPointers(const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ConstValues getValuesPointers(const size_type,
                                                          const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));

   private:
    //!
    std::array<ScalarComponentFunctionView, N> function_components;
  };

}  // namespace mgis::function

#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.ixx"

#endif /* LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX */
