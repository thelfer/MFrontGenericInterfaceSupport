/*!
 * \file   MGIS/Function/Tensors/CoalescedMemoryAccessCompositeTensorsView.hxx
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
#error "TFEL is required to use coalesced memory access tensor views"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX
#define LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX

#include <tuple>
#include <concepts>
#include <type_traits>
#include "TFEL/Math/Array/CoalescedView.hxx"
#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a tensor which is true if the type is equal to real or matches the
   * TensorConcept
   */
  template <typename T>
  concept ScalarOrTensorConcept = std::same_as<T, real> || TensorConcept<T>;

  /*!
   * \brief a coalescence view which acts as a tensorial function
   *
   * \tparam Space: functional space
   * \tparam N: number of components
   * \tparam is_mutable: boolean stating if the view can return mutable values.
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable = true>
  requires(N > 0)  //
      struct CoalescedMemoryAccessCompositeTensorsView
      : CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable> {
    //
    template <ScalarOrTensorConcept ValueType>
    using MutableValues =
        std::conditional_t<std::same_as<real, ValueType>,
                           real&,
                           tfel::math::CoalescedView<ValueType>>;
    //
    template <ScalarOrTensorConcept ValueType>
    using ConstValues =
        std::conditional_t<std::same_as<real, ValueType>,
                           const real&,
                           tfel::math::CoalescedView<const ValueType>>;
    // inheriting constructor
    using CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
        CoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type, const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(
        const size_type) const  //
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(const size_type,
                                                       const size_type) const
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));

  };  // end of CoalescedMemoryAccessCompositeTensorsView

}  // namespace mgis::function

#include "MGIS/Function/Tensors/CoalescedMemoryAccessCompositeTensorsView.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSCOMPOSITETENSORSVIEW_HXX \
        */
