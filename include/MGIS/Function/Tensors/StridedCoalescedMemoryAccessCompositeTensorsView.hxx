/*!
 * \file
 * MGIS/Function/Tensors/StridedCoalescedMemoryAccessCompositeTensorsView.hxx
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

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use strided coalesced memory access tensor views"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TENSORS_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX
#define LIB_MGIS_FUNCTION_TENSORS_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX

#include <tuple>
#include <concepts>
#include <type_traits>
#include "TFEL/Math/Array/StridedCoalescedView.hxx"
#include "MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  template <ScalarOrTensorConcept ValueType>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue {
    using type = tfel::math::StridedCoalescedView<ValueType>;
  };

  template <>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue<real> {
    using type = real&;
  };

  template <ScalarOrTensorConcept ValueType>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewConstValue {
    using type = tfel::math::StridedCoalescedView<const ValueType>;
  };

  template <>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewConstValue<real> {
    using type = const real&;
  };

  /*!
   * \brief a strided coalescence view which acts as a tensorial function
   *
   * \tparam Space: functional space
   * \tparam N: number of components
   * \tparam is_mutable: boolean stating if the view can return mutable values.
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable = true>
  requires(N > 0)  //
      struct StridedCoalescedMemoryAccessCompositeTensorsView
      : StridedCoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable> {
    //
    template <ScalarOrTensorConcept ValueType>
    using MutableValues =
        typename StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue<
            ValueType>::type;
    //
    template <ScalarOrTensorConcept ValueType>
    using ConstValues =
        typename StridedCoalescedMemoryAccessCompositeTensorsViewConstValue<
            ValueType>::type;
    // inheriting constructor
    using StridedCoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
        StridedCoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearElementSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type, const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearQuadratureSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(
        const size_type) const  //
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearElementSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(const size_type,
                                                       const size_type) const
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearQuadratureSpaceConcept<Space>);

  };  // end of StridedCoalescedMemoryAccessCompositeTensorsView

}  // namespace mgis::function

#include "MGIS/Function/Tensors/StridedCoalescedMemoryAccessCompositeTensorsView.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORS_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX \
        */
