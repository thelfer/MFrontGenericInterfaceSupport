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
      constexpr
      typename CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          MutableValues
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
          MutableValues
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
          ConstValues
      CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          getValuesPointers(const size_type i) const  //
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
          ConstValues
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

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_COALESCEDMEMORYACCESSFUNCTIONVIEWBASE_IXX */
