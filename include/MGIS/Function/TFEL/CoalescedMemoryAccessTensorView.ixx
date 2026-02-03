/*!
 * \file   MGIS/Function/TFEL/CoalescedMemoryAccessTensorView.ixx
 * \brief
 * \author Thomas Helfer
 * \date   27/10/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSTENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSTENSORVIEW_IXX

namespace mgis::function {

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename CoalescedMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::MutableValues
  CoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type
          i) requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                      (!hasElementWorkspace<Space>)) {
    return tfel::math::CoalescedView<TensorType>(this->getValuesPointers(i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename CoalescedMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::MutableValues
  CoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type e,
      const size_type q)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return tfel::math::CoalescedView<TensorType>(this->getValuesPointers(e, q));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename CoalescedMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::ConstValues
  CoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type i) const  //
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return tfel::math::CoalescedView<const TensorType>(
        this->getValuesPointers(i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename CoalescedMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::ConstValues
  CoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type e, const size_type q) const  //
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return tfel::math::CoalescedView<const TensorType>(
        this->getValuesPointers(e, q));
  }  // end of operator()

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSTENSORVIEW_IXX */
