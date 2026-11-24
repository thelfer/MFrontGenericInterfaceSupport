/*!
 * \file   MGIS/Function/StridedMemoryAccessFunctionViewBase.ixx
 * \brief File which contains the implementation of the
 * `StridedMemoryAccessFunctionViewBase` class declared in
 * `MGIS/Function/StridedMemoryAccessFunctionViewBase.hxx`.
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

#ifndef LIB_MGIS_FUNCTION_STRIDEDMEMORYACCESSFUNCTIONVIEWBASE_IXX
#define LIB_MGIS_FUNCTION_STRIDEDMEMORYACCESSFUNCTIONVIEWBASE_IXX

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr bool StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             std::span<const real> values)  //
      requires(!is_mutable) {
    const auto ne = getSpaceSize(s);
    if (values.size() != N * ne) {
      return eh.registerErrorMessage("invalid number of values");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          StridedMemoryAccessFunctionViewBase(
              const Space& s,
              std::span<const real> values) requires(!is_mutable)
      : StridedMemoryAccessFunctionViewBase(preconditions_check, s, values) {
  }  // end of StridedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          StridedMemoryAccessFunctionViewBase(
              const PreconditionsCheck<doPreconditionsCheck>& pcheck,
              const Space& s,
              std::span<const real> values) requires(!is_mutable)
      : PreconditionsChecker<StridedMemoryAccessFunctionViewBase>(
            pcheck, s, values),
        space(s),
        data_pointer(values.data()) {
  }  // end of StridedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr bool StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             std::span<real> values) {
    const auto ne = getSpaceSize(s);
    if (values.size() != N * ne) {
      return eh.registerErrorMessage("invalid number of values");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          StridedMemoryAccessFunctionViewBase(const Space& s,
                                              std::span<real> values)
      : StridedMemoryAccessFunctionViewBase(preconditions_check, s, values) {
  }  // end of StridedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
          StridedMemoryAccessFunctionViewBase(
              const PreconditionsCheck<doPreconditionsCheck>& pcheck,
              const Space& s,
              std::span<real> values)
      : PreconditionsChecker<StridedMemoryAccessFunctionViewBase>(
            pcheck, s, values),
        space(s),
        data_pointer(values.data()) {
  }  // end of StridedMemoryAccessFunctionViewBase

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      constexpr const Space& StridedMemoryAccessFunctionViewBase<
          Space,
          N,
          is_mutable>::getSpace() const MGIS_NOEXCEPT {
    return this->space;
  }  // end of getSpace

  // the values of the components are gathered using a pack of indices; the
  // offset of the value of the `k`-th component at a given integration point
  // is the sum of the offset of the integration point and of `k` times the
  // size of the space
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)                    //
      constexpr std::array<real, N>  //
      StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::getValues(
          const size_type i) const requires(LinearElementSpaceConcept<Space>) {
    return [ this, i ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real, N> {
      const auto ne = getSpaceSize(this->space);
      return {this->data_pointer[i + ne * Is]...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValues

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)                    //
      constexpr std::array<real, N>  //
      StridedMemoryAccessFunctionViewBase<Space, N, is_mutable>::getValues(
          const size_type e,
          const size_type q) const  //
      requires(LinearQuadratureSpaceConcept<Space>) {
    return [ this, e, q ]<std::size_t... Is>(std::index_sequence<Is...>)
        ->std::array<real, N> {
      const auto ne = getSpaceSize(this->space);
      const auto o = getQuadraturePointOffset(this->space, e, q);
      return {this->data_pointer[o + ne * Is]...};
    }
    (std::make_index_sequence<N>());
  }  // end of getValues

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_STRIDEDMEMORYACCESSFUNCTIONVIEWBASE_IXX */
