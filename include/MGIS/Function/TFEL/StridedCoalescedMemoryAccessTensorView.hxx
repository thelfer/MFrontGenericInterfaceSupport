/*!
 * \file   MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.hxx
 * \brief  File which defines the `StridedCoalescedMemoryAccessTensorView`
 * class, a view which acts as a function returning TFEL's tensorial objects
 * whose components are stored contiguously in a single memory region, in a
 * non interleaved manner.
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

    /*!
     * \brief partial specialization allowing functions to return a
     * `tfel::math::StridedCoalescedView`, i.e. a tensorial object whose
     * components are stored contiguously in a single memory region, in a non
     * interleaved manner
     *
     * \tparam TensorType: type of the tensorial object viewed, which must
     * match the `TensorConcept` concept
     */
    template <TensorConcept TensorType>
    struct FunctionResultTypeTraits<
        ::tfel::math::StridedCoalescedView<TensorType>> {
      //! \brief the traits are specialized for this type
      static constexpr auto is_specialized = true;
    };

  }  // namespace internals

  /*!
   * \brief a strided coalesced memory access view which acts as a function
   * returning TFEL's tensorial objects
   *
   * The view is built from a space and an array of values, whose number of
   * values must be equal to `N * getSpaceSize(space)`, `N` being the number of
   * components of `TensorType`. The values of each component are stored
   * contiguously in a single memory region, in a non interleaved manner, as
   * described by the `StridedCoalescedMemoryAccessFunctionViewBase` class: the
   * values of the `k`-th component start at the offset
   * `k * getSpaceSize(space)`, so the stride between the values of two
   * successive components at a given integration point is equal to the size of
   * the space.
   *
   * The call operators return a `tfel::math::StridedCoalescedView` of the
   * tensorial type, built from the address of the first component value at the
   * given integration point and the stride between successive components, and
   * which allows the tensor to be read, or written for mutable views, in a
   * single operation.
   *
   * Only the call operators which match the way the space can be addressed are
   * available.
   *
   * \tparam Space: functional space on which the view is defined
   * \tparam TensorType: type of the tensorial object mapped, which must match
   * the `TensorConcept` concept
   * \tparam is_mutable: boolean stating if the view can return mutable values
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable = true>
  struct StridedCoalescedMemoryAccessTensorView
      : StridedCoalescedMemoryAccessFunctionViewBase<
            Space,
            internals::CompileTimeSize<TensorType>::value,
            is_mutable> {
    //! \brief type returned by the mutable call operators
    using MutableValues = tfel::math::StridedCoalescedView<TensorType>;
    //! \brief type returned by the constant call operators
    using ConstValues = tfel::math::StridedCoalescedView<const TensorType>;

    /*!
     * \brief constructors inherited from the
     * `StridedCoalescedMemoryAccessFunctionViewBase` base class
     */
    using StridedCoalescedMemoryAccessFunctionViewBase<
        Space,
        internals::CompileTimeSize<TensorType>::value,
        is_mutable>::StridedCoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the values of all the components at the given integration point,
     * as a mutable strided coalesced view of the tensorial type
     * \param[in] i: offset associated with the integration point
     * \note only available for mutable views defined on a linear element space
     */
    [[nodiscard]] constexpr MutableValues operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space>);
    /*!
     * \return the values of all the components at the given integration point,
     * as a mutable strided coalesced view of the tensorial type
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for mutable views defined on a linear quadrature
     * space
     */
    [[nodiscard]] constexpr MutableValues
    operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space>);
    /*!
     * \return the values of all the components at the given integration point,
     * as a read-only strided coalesced view of the tensorial type
     * \param[in] i: offset associated with the integration point
     * \note only available for views defined on a linear element space
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space>);
    /*!
     * \return the values of all the components at the given integration point,
     * as a read-only strided coalesced view of the tensorial type
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available for views defined on a linear quadrature space
     */
    [[nodiscard]] constexpr ConstValues operator()(const size_type,
                                                   const size_type) const
        requires(LinearQuadratureSpaceConcept<Space>);

  };  // end of StridedCoalescedMemoryAccessTensorView

  /*!
   * \brief partial specialisation stating that a
   * `StridedCoalescedMemoryAccessTensorView` is its own view: it only stores
   * the space and a pointer to the first value, so copying it is cheap
   *
   * \tparam Space: functional space on which the view is defined
   * \tparam TensorType: type of the tensorial object mapped
   * \tparam is_mutable: boolean stating if the view can return mutable values
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  struct LightweightViewTraits<
      StridedCoalescedMemoryAccessTensorView<Space, TensorType, is_mutable>>
      : std::true_type {};

  /*!
   * \brief return the underlying space
   * \tparam Space: functional space on which the view is defined
   * \tparam TensorType: type of the tensorial object mapped
   * \tparam is_mutable: boolean stating if the view can return mutable values
   * \param[in] v: strided coalesced memory access tensor view
   * \return the space on which the view is defined
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr const Space& getSpace(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) MGIS_NOEXCEPT;
  /*!
   * \brief a noop function to match the `EvaluatorConcept` concept
   * \tparam Space: functional space on which the view is defined
   * \tparam TensorType: type of the tensorial object mapped
   * \tparam is_mutable: boolean stating if the view can return mutable values
   * \param[in] eh: unused error handler
   * \param[in] v: unused strided coalesced memory access tensor view
   * \return always true: the validity of the view is guaranteed by its
   * constructors
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) MGIS_NOEXCEPT;
  /*!
   * \brief return the number of components
   * \tparam Space: functional space on which the view is defined
   * \tparam TensorType: type of the tensorial object mapped
   * \tparam is_mutable: boolean stating if the view can return mutable values
   * \param[in] v: unused strided coalesced memory access tensor view
   * \return the number of components of the view, which is the number of
   * components of the tensorial type
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const StridedCoalescedMemoryAccessTensorView<Space,
                                                   TensorType,
                                                   is_mutable>&) MGIS_NOEXCEPT;

}  // namespace mgis::function

#include "MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSTENSORVIEW_HXX \
        */
