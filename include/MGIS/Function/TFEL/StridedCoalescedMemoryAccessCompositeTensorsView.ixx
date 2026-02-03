/*!
 * \file
 * MGIS/Function/TFEL/StridedCoalescedMemoryAccessCompositeTensorsView.ixx
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

#ifndef LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr typename StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::template MutableValues<ValueType>  //
      StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::get(const size_type i)  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               is_mutable && LinearElementSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    if constexpr (std::same_as<real, ValueType>) {
      return this->data_pointer[begin * ne + i];
    } else {
      return tfel::math::StridedCoalescedView<ValueType>(
          this->data_pointer + begin * ne + i, ne);
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr typename StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::template MutableValues<ValueType>  //
      StridedCoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          get(const size_type e,
              const size_type q)  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               is_mutable && LinearQuadratureSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    const auto o = getQuadraturePointOffset(this->space, e, q);
    if constexpr (std::same_as<real, ValueType>) {
      return this->data_pointer[begin * ne + o];
    } else {
      return tfel::math::StridedCoalescedView<ValueType>(
          this->data_pointer + begin * ne + o, ne);
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr typename StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::template ConstValues<ValueType>  //
      StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::get(const size_type i) const  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               LinearElementSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    if constexpr (std::same_as<real, ValueType>) {
      return this->data_pointer[begin * ne + i];
    } else {
      return tfel::math::StridedCoalescedView<const ValueType>(
          this->data_pointer + begin * ne + i, ne);
    }
  }  // end of get

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)  //
      template <size_type begin, ScalarOrTensorConcept ValueType>
      constexpr typename StridedCoalescedMemoryAccessCompositeTensorsView<
          Space,
          N,
          is_mutable>::template ConstValues<ValueType>  //
      StridedCoalescedMemoryAccessCompositeTensorsView<Space, N, is_mutable>::
          get(const size_type e,
              const size_type q) const  //
      requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
               LinearQuadratureSpaceConcept<Space>) {
    const auto ne = getSpaceSize(this->space);
    const auto o = getQuadraturePointOffset(this->space, e, q);
    if constexpr (std::same_as<real, ValueType>) {
      return this->data_pointer[begin * ne + o];
    } else {
      return tfel::math::StridedCoalescedView<const ValueType>(
          this->data_pointer + begin * ne + o, ne);
    }
  }  // end of get

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_IXX \
        */
