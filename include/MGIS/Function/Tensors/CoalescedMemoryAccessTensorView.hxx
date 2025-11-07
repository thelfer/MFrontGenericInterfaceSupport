/*!
 * \file   MGIS/Function/Tensors/CoalescedMemoryAccessTensorView.hxx
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

#ifndef LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSTENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSTENSORVIEW_HXX

#include "TFEL/Math/Array/CoalescedView.hxx"
#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a coalescence view which acts as a tensorial function
   *
   * \tparam Space: functional space
   * \tparam TensorType: tensorial object mapped
   * \tparam is_mutable: boolean stating if the view can return mutable values.
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable = true>
  struct CoalescedMemoryAccessTensorView
      : CoalescedMemoryAccessFunctionViewBase<
            Space,
            internals::CompileTimeSize<TensorType>::value,
            is_mutable> {
    //
    using MutableValues = tfel::math::CoalescedView<TensorType>;
    //
    using ConstValues = tfel::math::CoalescedView<const TensorType>;

    // inheriting constructor
    using CoalescedMemoryAccessFunctionViewBase<
        Space,
        internals::CompileTimeSize<TensorType>::value,
        is_mutable>::CoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr MutableValues operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr MutableValues
    operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type,
                                                   const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));

  };  // end of CoalescedMemoryAccessTensorView

}  // namespace mgis::function

#include "MGIS/Function/Tensors/CoalescedMemoryAccessTensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORS_COALESCEDMEMORYACCESSTENSORVIEW_HXX */
