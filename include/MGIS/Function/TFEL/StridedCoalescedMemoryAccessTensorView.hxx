/*!
 * \file   MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.hxx
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

#ifndef LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_HXX

#include "TFEL/Math/Array/StridedCoalescedView.hxx"
#include "MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

namespace mgis::function {

  namespace internals {

    //! \brief partial specialization for strided coalesced views
    template <TensorConcept TensorType>
    struct FunctionResultTypeTraits<
        ::tfel::math::StridedCoalescedView<TensorType>> {
      static constexpr auto is_specialized = true;
    };

  }  // namespace internals

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
  struct StridedCoalescedMemoryAccessTensorView
      : StridedCoalescedMemoryAccessFunctionViewBase<
            Space,
            internals::CompileTimeSize<TensorType>::value,
            is_mutable> {
    //
    using MutableValues = tfel::math::StridedCoalescedView<TensorType>;
    //
    using ConstValues = tfel::math::StridedCoalescedView<const TensorType>;

    // inheriting constructor
    using StridedCoalescedMemoryAccessFunctionViewBase<
        Space,
        internals::CompileTimeSize<TensorType>::value,
        is_mutable>::StridedCoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr MutableValues operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr MutableValues
    operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space>);
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type,
                                                   const size_type) const
        requires(LinearQuadratureSpaceConcept<Space>);

  };  // end of StridedCoalescedMemoryAccessTensorView

  //! \brief partial specialisation
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  struct LightweightViewTraits<
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>>
      : std::true_type {};

  //! \return the underlying space
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr const Space& getSpace(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) noexcept;
  //! \brief a noop function to match the EvaluatorConcept concept
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) noexcept;
  //! \return the number of components
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) noexcept;

}  // namespace mgis::function

#include "MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_HXX \
        */
