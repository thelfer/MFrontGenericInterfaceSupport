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

  template <size_type N, typename Space>
  constexpr std::optional<
      std::array<FunctionView<Space,
                              FunctionDataLayoutDescription{.data_size = 1,
                                                            .data_stride = 1}>,
                 N>>
  splitArrayIntoScalarFunctionViews(AbstractErrorHandler&,
                                    const Space&,
                                    std::span<real>);

  template <size_type N, typename Space>
  constexpr std::optional<
      std::array<FunctionView<Space,
                              FunctionDataLayoutDescription{.data_size = 1,
                                                            .data_stride = 1},
                              false>,
                 N>>
  splitArrayIntoScalarFunctionViews(AbstractErrorHandler&,
                                    const Space&,
                                    std::span<const real>);

  /*!
   *
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) struct CoalescedMemoryAccessFunctionViewBase
      : private PreconditionsChecker<
            CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>> {
    //
    template <size_type size>
    using MutableValues = std::array<real*, size>;
    //
    template <size_type size>
    using ConstValues = std::array<const real*, size>;

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
    constexpr CoalescedMemoryAccessFunctionViewBase(
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
    constexpr CoalescedMemoryAccessFunctionViewBase(const Space&,
                                                    std::span<real>);
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr MutableValues<N> getValuesPointers(
        const size_type)  //
        requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr MutableValues<N> getValuesPointers(
        const size_type, const size_type)  //
        requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated wNith the integration point
     */
    [[nodiscard]] constexpr ConstValues<N> getValuesPointers(
        const size_type) const  //
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ConstValues<N> getValuesPointers(
        const size_type, const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr MutableValues<size> getValuesPointers(
        const size_type)  //
        requires((begin + size <= N) && is_mutable &&
                 LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr MutableValues<size> getValuesPointers(
        const size_type,
        const size_type)  //
        requires((begin + size <= N) && is_mutable &&
                 LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr ConstValues<size> getValuesPointers(
        const size_type) const
        requires((begin + size <= N) && LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr ConstValues<size> getValuesPointers(
        const size_type, const size_type) const
        requires((begin + size <= N) && LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type offset>
    [[nodiscard]] constexpr real* getValuePointer(const size_type)  //
        requires((offset < N) && is_mutable &&
                 LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type offset>
    [[nodiscard]] constexpr real* getValuePointer(const size_type,
                                                  const size_type)  //
        requires((offset < N) && is_mutable &&
                 LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type offset>
    [[nodiscard]] constexpr const real* getValuePointer(const size_type) const
        requires((offset < N) && LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type offset>
    [[nodiscard]] constexpr const real* getValuePointer(const size_type,
                                                        const size_type) const
        requires((offset < N) && LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));

   protected:
    //
    static constexpr auto splitValues(const Space&, std::span<real>);
    //
    static constexpr auto splitValues(const Space&, std::span<const real>);
    /*!
     * \param[in] components: components
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \param[in] values: values
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<const real>) requires(!is_mutable);
    /*!
     * \param[in] values: values
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const Space&,
        std::span<real>);
    //! \brief function's view for each component
    std::array<ScalarComponentFunctionView, N> function_components;
  };

}  // namespace mgis::function

#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.ixx"

#endif /* LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_HXX */
