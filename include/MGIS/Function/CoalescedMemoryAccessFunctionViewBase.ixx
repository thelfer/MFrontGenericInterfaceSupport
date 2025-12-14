/*!
 * \file   MGIS/Function/CoalescedMemoryAccessFunctionViewBase.ixx
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

#ifndef LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_IXX
#define LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_IXX

namespace mgis::function {

  template <size_type N, typename Space>
  constexpr std::optional<
      std::array<FunctionView<Space,
                              FunctionDataLayoutDescription{.data_size = 1,
                                                            .data_stride = 1}>,
                 N>>
  splitArrayIntoScalarFunctionViews(AbstractErrorHandler& ctx,
                                    const Space& space,
                                    std::span<real> values) {
    using ScalarFunctionView =
        FunctionView<Space, FunctionDataLayoutDescription{.data_size = 1,
                                                          .data_stride = 1}>;
    const auto ne = getSpaceSize(space);
    if (values.size() != N * ne) {
      return ctx.registerErrorMessage("invalid number of values");
    }
    return
        [&space, &values, &ne ]<std::size_t... Is>(std::index_sequence<Is...>)
            ->std::array<ScalarFunctionView, N> {
      return {ScalarFunctionView(space, values.subspan(Is * ne, ne))...};
    }
    (std::make_index_sequence<N>());
  }  // end of splitArrayIntoScalarFunctionViews

  template <size_type N, typename Space>
  constexpr std::optional<
      std::array<FunctionView<Space,
                              FunctionDataLayoutDescription{.data_size = 1,
                                                            .data_stride = 1},
                              false>,
                 N>>
  splitArrayIntoScalarFunctionViews(AbstractErrorHandler& ctx,
                                    const Space& space,
                                    std::span<const real> values) {
    using ScalarFunctionView = FunctionView<
        Space, FunctionDataLayoutDescription{.data_size = 1, .data_stride = 1},
        false>;
    const auto ne = getSpaceSize(space);
    if (values.size() != N * ne) {
      return ctx.registerErrorMessage("invalid number of values");
    }
    return
        [&space, &values, &ne ]<std::size_t... Is>(std::index_sequence<Is...>)
            ->std::array<ScalarFunctionView, N> {
      return {ScalarFunctionView(space, values.subspan(Is * ne, ne))...};
    }
    (std::make_index_sequence<N>());
  }  // end of splitArrayIntoScalarFunctionViews

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr auto CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::splitValues(const Space& space, std::span<real> values) {
    auto ctx = ContractViolationHandler{};
    return *(splitArrayIntoScalarFunctionViews<N>(ctx, space, values));
  }  // end of splitValues

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr auto CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::splitValues(const Space& space,
                                   std::span<const real> values) {
    auto ctx = ContractViolationHandler{};
    return *(splitArrayIntoScalarFunctionViews<N>(ctx, space, values));
  }  // end of splitValues

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr bool CoalescedMemoryAccessFunctionViewBase<Space,
                                                           N,
                                                           is_mutable>::
          checkPreconditions(
              AbstractErrorHandler& eh,
              const std::array<ScalarComponentFunctionView, N>& components) {
    for (int i = 1; i != N; ++i) {
      if (!areEquivalent(getSpace(components[0]), getSpace(components[1]))) {
        return eh.registerErrorMessage("unmatched spaces");
      }
    }
    for (const auto& c : components) {
      if (!check(eh, c)) {
        return false;
      }
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(
              const std::array<ScalarComponentFunctionView, N>& components)
      : CoalescedMemoryAccessFunctionViewBase(preconditions_check, components) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(
              const PreconditionsCheck<doPreconditionsCheck>& pcheck,
              const std::array<ScalarComponentFunctionView, N>& components)
      : PreconditionsChecker<CoalescedMemoryAccessFunctionViewBase>(pcheck,
                                                                    components),
        function_components(components) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr bool CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::checkPreconditions(AbstractErrorHandler& eh,
                                          const Space& space,
                                          std::span<const real> values)  //
      requires(!is_mutable) {
    const auto ne = getSpaceSize(space);
    if (values.size() != N * ne) {
      return eh.registerErrorMessage("invalid number of values");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(
              const Space& space,
              std::span<const real> values) requires(!is_mutable)
      : CoalescedMemoryAccessFunctionViewBase(
            preconditions_check, space, values) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(
              const PreconditionsCheck<doPreconditionsCheck>& pcheck,
              const Space& space,
              std::span<const real> values) requires(!is_mutable)
      : PreconditionsChecker<CoalescedMemoryAccessFunctionViewBase>(
            pcheck, space, values),
        function_components(splitValues(space, values)) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr bool CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::checkPreconditions(AbstractErrorHandler& eh,
                                          const Space& space,
                                          std::span<real> values) {
    const auto ne = getSpaceSize(space);
    if (values.size() != N * ne) {
      return eh.registerErrorMessage("invalid number of values");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(const Space& space,
                                                std::span<real> values)
      : CoalescedMemoryAccessFunctionViewBase(
            preconditions_check, space, values) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          CoalescedMemoryAccessFunctionViewBase(
              const PreconditionsCheck<doPreconditionsCheck>& pcheck,
              const Space& space,
              std::span<real> values)
      : PreconditionsChecker<CoalescedMemoryAccessFunctionViewBase>(
            pcheck, space, values),
        function_components(splitValues(space, values)) {
  }  // end of CoalescedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          MutableValues<>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type i) requires(
              is_mutable&& LinearElementSpaceConcept<Space> &&
              (!hasElementWorkspace<Space>)) {
    return [ this, i ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, N> {
      return {&this->function_components[Is](i)...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          MutableValues<>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type e,
                            const size_type q)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return [ this, e, q ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, N> {
      return {&this->function_components[Is](e, q)...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          ConstValues<>  //
      CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::getValuesPointers(const size_type i) const  //
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return [ this, i ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<const real*, N> {
      return {&this->function_components[Is](i)...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          ConstValues<>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type e, const size_type q) const  //
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return [ this, e, q ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, N> {
      return {&this->function_components[Is](e, q)...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, size_type size>
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          MutableValues<size>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type i)  //
      requires((begin + size <= N) && is_mutable &&
               LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return [ this, i ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, size> {
      return {&this->function_components[begin + Is](i)...};
    }
    (std::make_index_sequence<size>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, size_type size>
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          MutableValues<size>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type e,
                            const size_type q)  //
      requires((begin + size <= N) && is_mutable &&
               LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return [ this, e, q ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, size> {
      return {&this->function_components[begin + Is](e, q)...};
    }
    (std::make_index_sequence<size>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, size_type size>
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          ConstValues<size>  //
      CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::getValuesPointers(const size_type i) const  //
      requires((begin + size <= N) && LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return [ this, i ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<const real*, N> {
      return {&this->function_components[begin + Is](i)...};
    }
    (std::make_index_sequence<size>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, size_type size>
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          ConstValues<size>  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type e, const size_type q) const  //
      requires((begin + size <= N) && LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return [ this, e, q ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real*, size> {
      return {&this->function_components[begin + Is](e, q)...};
    }
    (std::make_index_sequence<size>());
  }  // end of getValuesPointers

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) template <size_type offset>
  constexpr real* CoalescedMemoryAccessFunctionViewBase<
      Space,
      N,
      is_mutable>::getValuePointer(const size_type i)  //
      requires((offset < N) && is_mutable && LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return &(this->function_components[offset](i));
  }  // end of getValuePointer

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type offset>
      constexpr real* CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::getValuePointer(const size_type e,
                                       const size_type q)  //
      requires((offset < N) && is_mutable &&
               LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return &(this->function_components[offset](e, q));
  }  // end of getValuePointer

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type offset>
      constexpr const real* CoalescedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::getValuePointer(const size_type i) const  //
      requires((offset < N) && LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return &(this->function_components[offset](i));
  }  // end of getValuePointer

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) template <size_type offset>
  constexpr const real*  //
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuePointer(const size_type e, const size_type q) const  //
      requires((offset < N) && LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return &(this->function_components[offset](e, q));
  }  // end of getValuePointer

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_IXX */
