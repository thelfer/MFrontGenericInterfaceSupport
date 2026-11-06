/*!
 * \file   MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.ixx
 * \brief  File which contains the implementation of the
 * `StridedCoalescedMemoryAccessTensorView` class declared in
 * `MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.hxx`.
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

#ifndef LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX

namespace mgis::function {

  // the call operators build a strided coalesced view from the address of the
  // first component value at the given integration point, the stride between
  // successive components being equal to the size of the space
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr
      typename StridedCoalescedMemoryAccessTensorView<Space,
                                                      TensorType,
                                                      is_mutable>::MutableValues
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::
      operator()(const size_type i) requires(
          is_mutable&& LinearElementSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    return tfel::math::StridedCoalescedView<TensorType>(this->data_pointer + i,
                                                        ne);
  }  // end of operator()

  // for a quadrature space, the offset of the quadrature point is first
  // computed using the `getQuadraturePointOffset` function
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr
      typename StridedCoalescedMemoryAccessTensorView<Space,
                                                      TensorType,
                                                      is_mutable>::MutableValues
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::
      operator()(const size_type e,
                 const size_type q)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    const auto o = getQuadraturePointOffset(this->space, e, q);
    return tfel::math::StridedCoalescedView<TensorType>(this->data_pointer + o,
                                                        ne);
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr
      typename StridedCoalescedMemoryAccessTensorView<Space,
                                                      TensorType,
                                                      is_mutable>::ConstValues
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::
      operator()(const size_type i) const  //
      requires(LinearElementSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    return tfel::math::StridedCoalescedView<const TensorType>(
        this->data_pointer + i, ne);
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr
      typename StridedCoalescedMemoryAccessTensorView<Space,
                                                      TensorType,
                                                      is_mutable>::ConstValues
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::
      operator()(const size_type e, const size_type q) const  //
      requires(LinearQuadratureSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    const auto o = getQuadraturePointOffset(this->space, e, q);
    return tfel::math::StridedCoalescedView<const TensorType>(
        this->data_pointer + o, ne);
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr const Space& getSpace(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>& v) MGIS_NOEXCEPT {
    return v.getSpace();
  }  // end of getSpace

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr bool check(
      AbstractErrorHandler&,
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) MGIS_NOEXCEPT {
    return true;
  }  // end of check

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr size_type getNumberOfComponents(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) MGIS_NOEXCEPT {
    return compile_time_size<TensorType>;
  }  // end of getNumberOfComponents

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX \
        */
