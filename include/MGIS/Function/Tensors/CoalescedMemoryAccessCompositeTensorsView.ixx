/*!
 * \file   MGIS/Function/Tensors/CoalescedMemoryAccessCompositeTensorsView.ixx
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

#ifndef LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX
#define LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr
      typename CoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          template MutableValues<ValueType> CoalescedMemoryAccessCompositeTensorsView<
              Space,
              N,
              is_mutable>::get(const size_type i)  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               is_mutable && LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    if constexpr (std::same_as<real, ValueType>) {
      const auto ptr = this->template getValuePointer<begin>(i);
      return *(ptr[0]);
    } else {
      return tfel::math::CoalescedView<ValueType>(
          this->template getValuesPointers<
              begin, internals::CompileTimeSize<ValueType>::value>(i));
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr
      typename CoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          template MutableValues<ValueType> CoalescedMemoryAccessCompositeTensorsView<
              Space,
              N,
              is_mutable>::get(const size_type e,
                               const size_type q)  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               is_mutable && LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    if constexpr (std::same_as<real, ValueType>) {
      const auto ptr = this->template getValuePointer<begin>(e, q);
      return *(ptr[0]);
    } else {
      return tfel::math::CoalescedView<ValueType>(
          this->template getValuesPointers<
              begin, internals::CompileTimeSize<ValueType>::value>(e, q));
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr
      typename CoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          template ConstValues<ValueType> CoalescedMemoryAccessCompositeTensorsView<
              Space,
              N,
              is_mutable>::get(const size_type i) const  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    if constexpr (std::same_as<real, ValueType>) {
      const auto ptr = this->template getValuePointer<begin>(i);
      return *(ptr[0]);
    } else {
      return tfel::math::CoalescedView<const ValueType>(
          this->template getValuesPointers<
              begin, internals::CompileTimeSize<ValueType>::value>(i));
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr
      typename CoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          template ConstValues<ValueType> CoalescedMemoryAccessCompositeTensorsView<
              Space,
              N,
              is_mutable>::get(const size_type e,
                               const size_type q) const  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    if constexpr (std::same_as<real, ValueType>) {
      const auto ptr = this->template getValuePointer<begin>(e, q);
      return *(ptr[0]);
    } else {
      return tfel::math::CoalescedView<const ValueType>(
          this->template getValuesPointers<
              begin, internals::CompileTimeSize<ValueType>::value>(e, q));
    }
  }  // end of get

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX \
        */
