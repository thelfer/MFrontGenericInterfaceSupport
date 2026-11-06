/*!
 * \file   MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx
 * \brief File which defines a base class for views whose components are
 * stored contiguously in separate memory regions, each component being
 * described by a scalar function view.
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

  /*!
   * \brief split an array of values into `N` mutable scalar function views,
   * each of them being defined on the given space
   *
   * The values of the `k`-th component are stored in
   * `values.subspan(k * getSpaceSize(space), getSpaceSize(space))`, i.e. each
   * component is stored contiguously and the components are not interleaved.
   *
   * \tparam N: number of components
   * \tparam Space: type of the space on which the views are defined
   * \param[in] ctx: error handler used to report the errors
   * \param[in] space: space on which the views are defined
   * \param[in] values: values
   * \return the views of the components or an empty optional if the number of
   * values is not equal to `N * getSpaceSize(space)`
   */
  template <size_type N, typename Space>
  constexpr std::optional<
      std::array<FunctionView<Space,
                              FunctionDataLayoutDescription{.data_size = 1,
                                                            .data_stride = 1}>,
                 N>>
  splitArrayIntoScalarFunctionViews(AbstractErrorHandler&,
                                    const Space&,
                                    std::span<real>);
  /*!
   * \brief split an array of values into `N` constant scalar function views,
   * each of them being defined on the given space
   *
   * The values of the `k`-th component are stored in
   * `values.subspan(k * getSpaceSize(space), getSpaceSize(space))`, i.e. each
   * component is stored contiguously and the components are not interleaved.
   *
   * \tparam N: number of components
   * \tparam Space: type of the space on which the views are defined
   * \param[in] ctx: error handler used to report the errors
   * \param[in] space: space on which the views are defined
   * \param[in] values: constant values
   * \return the constant views of the components or an empty optional if the
   * number of values is not equal to `N * getSpaceSize(space)`
   */
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
   * \brief base class for views whose components are stored contiguously in
   * separate memory regions
   *
   * Each component is described by a scalar function view whose data are
   * stored contiguously for the different points of the space, as illustrated
   * by the following scheme:
   *
   * ~~~~
   * | <------- Component 1 ---------> |....| <----- Component N ---------> |
   * +-------++-------++------++-------+----+-------++------++------++-------+
   * | Elt 1 || Elt 2 || .... || Elt N |....| Elt 1 ||Elt 2 || .... || Elt N |
   * +-------++-------++------++-------+----+-------++------++------++-------+
   * ~~~~
   *
   * The values of the different components at a given integration point are
   * not contiguous: the `getValuesPointers` function gathers the addresses of
   * those values, which allows the values of all the components at a given
   * integration point to be read or written in a single operation. This
   * scheme is used by the `CoalescedMemoryAccessTensorView` class to map
   * tensorial objects whose components are stored in separate contiguous
   * memory regions.
   *
   * \tparam Space: type of the space on which the views are defined
   * \tparam N: number of components, which must be strictly positive
   * \tparam is_mutable: boolean stating if the view can return mutable pointers
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) struct CoalescedMemoryAccessFunctionViewBase
      : private PreconditionsChecker<
            CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>> {
    /*!
     * \brief array of the addresses of the values of `size` mutable components
     * at a given integration point
     * \tparam size: number of components
     */
    template <size_type size>
    using MutableValues = std::array<real*, size>;
    /*!
     * \brief array of the addresses of the values of `size` constant
     * components at a given integration point
     * \tparam size: number of components
     */
    template <size_type size>
    using ConstValues = std::array<const real*, size>;

    /*!
     * \brief type of the function view associated with a single component,
     * whose values are stored contiguously for the different points of the
     * space
     */
    using ScalarComponentFunctionView =
        FunctionView<Space,
                     FunctionDataLayoutDescription{.data_size = 1,
                                                   .data_stride = 1},
                     is_mutable>;
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] components: views of the components
     * \return true if the components are defined on equivalent spaces and if
     * all of them are valid
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \brief constructor, checking that the components share the same space
     * \param[in] components: views of the components
     */
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const std::array<ScalarComponentFunctionView, N>&);
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
     * \brief constructor splitting the values into `N` constant scalar
     * function views, available only for constant views
     * \param[in] space: space on which the view is defined
     * \param[in] values: constant values
     */
    constexpr CoalescedMemoryAccessFunctionViewBase(
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
     * \brief constructor splitting the values into `N` scalar function views
     * \param[in] space: space on which the view is defined
     * \param[in] values: values
     */
    constexpr CoalescedMemoryAccessFunctionViewBase(const Space&,
                                                    std::span<real>);
    /*!
     * \brief return the underlying space
     * \return the space on which the view is defined, i.e. the space of the
     * first component
     */
    [[nodiscard]] constexpr const Space& getSpace() const MGIS_NOEXCEPT;
    /*!
     * \return the addresses of the values of all the components at the given
     * integration point
     * \param[in] i: offset associated with the integration point
     * \note only available for mutable views defined on a linear element space
     * without element workspace
     */
    [[nodiscard]] constexpr MutableValues<N> getValuesPointers(
        const size_type)  //
        requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the addresses of the values of all the components at the given
     * integration point
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for mutable views defined on a linear quadrature
     * space without cell workspace
     */
    [[nodiscard]] constexpr MutableValues<N> getValuesPointers(
        const size_type, const size_type)  //
        requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the addresses of the values of all the components at the given
     * integration point
     * \param[in] i: offset associated with the integration point
     * \note only available for views defined on a linear element space without
     * element workspace
     */
    [[nodiscard]] constexpr ConstValues<N> getValuesPointers(
        const size_type) const  //
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the addresses of the values of all the components at the given
     * integration point
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for views defined on a linear quadrature space
     * without cell workspace
     */
    [[nodiscard]] constexpr ConstValues<N> getValuesPointers(
        const size_type, const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the addresses of the values of the components with indices from
     * `begin` to `begin + size - 1` at the given integration point
     * \tparam begin: index of the first component
     * \tparam size: number of components, such that `begin + size` is not
     * greater than `N`
     * \param[in] i: offset associated with the integration point
     * \note only available for mutable views defined on a linear element space
     * without element workspace
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr MutableValues<size> getValuesPointers(
        const size_type)  //
        requires((begin + size <= N) && is_mutable &&
                 LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the addresses of the values of the components with indices from
     * `begin` to `begin + size - 1` at the given integration point
     * \tparam begin: index of the first component
     * \tparam size: number of components, such that `begin + size` is not
     * greater than `N`
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for mutable views defined on a linear quadrature
     * space without cell workspace
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr MutableValues<size> getValuesPointers(
        const size_type,
        const size_type)  //
        requires((begin + size <= N) && is_mutable &&
                 LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the addresses of the values of the components with indices from
     * `begin` to `begin + size - 1` at the given integration point
     * \tparam begin: index of the first component
     * \tparam size: number of components, such that `begin + size` is not
     * greater than `N`
     * \param[in] i: offset associated with the integration point
     * \note only available for views defined on a linear element space without
     * element workspace
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr ConstValues<size> getValuesPointers(
        const size_type) const
        requires((begin + size <= N) && LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the addresses of the values of the components with indices from
     * `begin` to `begin + size - 1` at the given integration point
     * \tparam begin: index of the first component
     * \tparam size: number of components, such that `begin + size` is not
     * greater than `N`
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for views defined on a linear quadrature space
     * without cell workspace
     */
    template <size_type begin, size_type size>
    [[nodiscard]] constexpr ConstValues<size> getValuesPointers(
        const size_type, const size_type) const
        requires((begin + size <= N) && LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the address of the value of the `offset`-th component at the
     * given integration point
     * \tparam offset: index of the component, strictly smaller than `N`
     * \param[in] i: offset associated with the integration point
     * \note only available for mutable views defined on a linear element space
     * without element workspace
     */
    template <size_type offset>
    [[nodiscard]] constexpr real* getValuePointer(const size_type)  //
        requires((offset < N) && is_mutable &&
                 LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the address of the value of the `offset`-th component at the
     * given integration point
     * \tparam offset: index of the component, strictly smaller than `N`
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for mutable views defined on a linear quadrature
     * space without cell workspace
     */
    template <size_type offset>
    [[nodiscard]] constexpr real* getValuePointer(const size_type,
                                                  const size_type)  //
        requires((offset < N) && is_mutable &&
                 LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the address of the value of the `offset`-th component at the
     * given integration point
     * \tparam offset: index of the component, strictly smaller than `N`
     * \param[in] i: offset associated with the integration point
     * \note only available for views defined on a linear element space without
     * element workspace
     */
    template <size_type offset>
    [[nodiscard]] constexpr const real* getValuePointer(const size_type) const
        requires((offset < N) && LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the address of the value of the `offset`-th component at the
     * given integration point
     * \tparam offset: index of the component, strictly smaller than `N`
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for views defined on a linear quadrature space
     * without cell workspace
     */
    template <size_type offset>
    [[nodiscard]] constexpr const real* getValuePointer(const size_type,
                                                        const size_type) const
        requires((offset < N) && LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));

   protected:
    /*!
     * \brief split the values into `N` mutable scalar function views
     * \param[in] space: space on which the views are defined
     * \param[in] values: values, whose number of values must be equal to
     * `N * getSpaceSize(space)`
     * \return the views of the components
     * \note this function assumes that the preconditions have been checked
     * before being called
     */
    static constexpr auto splitValues(const Space&, std::span<real>);
    /*!
     * \brief split the values into `N` constant scalar function views
     * \param[in] space: space on which the views are defined
     * \param[in] values: constant values, whose number of values must be
     * equal to `N * getSpaceSize(space)`
     * \return the constant views of the components
     * \note this function assumes that the preconditions have been checked
     * before being called
     */
    static constexpr auto splitValues(const Space&, std::span<const real>);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] components: views of the components
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const std::array<ScalarComponentFunctionView, N>&);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] space: space on which the view is defined
     * \param[in] values: constant values
     */
    template <bool doPreconditionsCheck>
    constexpr CoalescedMemoryAccessFunctionViewBase(
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
