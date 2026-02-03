/*!
 * \file   MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.ixx
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

#ifndef LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX

namespace mgis::function {

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

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_IXX \
        */
