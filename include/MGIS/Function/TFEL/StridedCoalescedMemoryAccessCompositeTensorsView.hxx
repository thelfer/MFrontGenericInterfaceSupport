/*!
 * \file MGIS/Function/TFEL/StridedCoalescedMemoryAccessCompositeTensorsView.hxx
 * \brief  File which defines the
 * `StridedCoalescedMemoryAccessCompositeTensorsView` class, a view handling
 * several tensorial objects of possibly different types sharing the same set
 * of components, whose values are stored contiguously in a single memory
 * region, in a non interleaved manner.
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

#ifndef LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX

#include <tuple>
#include <concepts>
#include <type_traits>
#include "TFEL/Math/Array/StridedCoalescedView.hxx"
#include "MGIS/Function/StridedCoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief traits defining the type used by the
   * `StridedCoalescedMemoryAccessCompositeTensorsView` class to return a
   * mutable value of the given type
   *
   * The primary template returns a strided coalesced view of the tensorial
   * type, whose components are stored contiguously in a single memory region,
   * in a non interleaved manner. It is specialized below for scalar values,
   * which are returned by reference.
   *
   * \tparam ValueType: type of the value to be returned, which must match the
   * `ScalarOrTensorConcept` concept
   */
  template <ScalarOrTensorConcept ValueType>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue {
    //! \brief type of the mutable value returned
    using type = tfel::math::StridedCoalescedView<ValueType>;
  };

  /*!
   * \brief specialization for scalar values, which are returned as mutable
   * references to the value of the corresponding component
   */
  template <>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue<real> {
    //! \brief type of the mutable value returned
    using type = real&;
  };

  /*!
   * \brief traits defining the type used by the
   * `StridedCoalescedMemoryAccessCompositeTensorsView` class to return a
   * read-only value of the given type
   *
   * The primary template returns a strided coalesced view of the const
   * tensorial type, whose components are stored contiguously in a single memory
   * region, in a non interleaved manner. It is specialized below for scalar
   * values, which are returned by const reference.
   *
   * \tparam ValueType: type of the value to be returned, which must match the
   * `ScalarOrTensorConcept` concept
   */
  template <ScalarOrTensorConcept ValueType>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewConstValue {
    //! \brief type of the read-only value returned
    using type = tfel::math::StridedCoalescedView<const ValueType>;
  };

  /*!
   * \brief specialization for scalar values, which are returned as read-only
   * references to the value of the corresponding component
   */
  template <>
  struct StridedCoalescedMemoryAccessCompositeTensorsViewConstValue<real> {
    //! \brief type of the read-only value returned
    using type = const real&;
  };

  /*!
   * \brief a strided coalesced memory access view handling several tensorial
   * objects of possibly different types sharing the same set of components
   *
   * The view is built from a space and an array of values, whose number of
   * values must be equal to `N * getSpaceSize(space)`, as described by the
   * `StridedCoalescedMemoryAccessFunctionViewBase` class. The values of each
   * component are stored contiguously in a single memory region, in a non
   * interleaved manner: the values of the `k`-th component start at the offset
   * `k * getSpaceSize(space)`.
   *
   * The `get` member function returns, for a given integration point, the
   * tensorial object built from the components with indices from `begin` to
   * `begin + M - 1`, `M` being the number of components of the queried type,
   * as a strided coalesced view whose stride is equal to the size of the space.
   * Several tensorial objects and scalar values may thus be extracted from the
   * same view, provided that their ranges of components do not exceed the `N`
   * available components, as stated by the `requires`-clause.
   *
   * Only the member functions which match the way the space can be addressed
   * are available.
   *
   * \tparam Space: functional space on which the view is defined
   * \tparam N: number of components, which must be strictly positive
   * \tparam is_mutable: boolean stating if the view can return mutable values
   */
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable = true>
  requires(N > 0)  //
      struct StridedCoalescedMemoryAccessCompositeTensorsView
      : StridedCoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable> {
    //
    /*!
     * \brief type returned by the mutable `get` member functions
     *
     * \tparam ValueType: type of the tensorial object or scalar value returned
     */
    template <ScalarOrTensorConcept ValueType>
    using MutableValues =
        typename StridedCoalescedMemoryAccessCompositeTensorsViewMutableValue<
            ValueType>::type;
    //
    /*!
     * \brief type returned by the constant `get` member functions
     *
     * \tparam ValueType: type of the tensorial object or scalar value returned
     */
    template <ScalarOrTensorConcept ValueType>
    using ConstValues =
        typename StridedCoalescedMemoryAccessCompositeTensorsViewConstValue<
            ValueType>::type;
    /*!
     * \brief constructors inherited from the
     * `StridedCoalescedMemoryAccessFunctionViewBase` base class
     */
    using StridedCoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
        StridedCoalescedMemoryAccessFunctionViewBase;
    /*!
     * \return the tensorial object built from the components of the view at the
     * given integration point
     * \tparam begin: index of the first component of the tensorial object
     * \tparam ValueType: type of the tensorial object or scalar value returned
     * \param[in] i: offset associated with the integration point
     * \note only available if the components with indices from `begin` to
     * `begin + M - 1` exist, `M` being the number of components of
     * `ValueType`, and for mutable views defined on a linear element space
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearElementSpaceConcept<Space>);
    /*!
     * \return the tensorial object built from the components of the view at the
     * given integration point
     * \tparam begin: index of the first component of the tensorial object
     * \tparam ValueType: type of the tensorial object or scalar value returned
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available if the components with indices from `begin` to
     * `begin + M - 1` exist, `M` being the number of components of
     * `ValueType`, and for mutable views defined on a linear quadrature space
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr MutableValues<ValueType>
    get(const size_type, const size_type) requires(
        (begin + internals::CompileTimeSize<ValueType>::value <= N) &&
        is_mutable && LinearQuadratureSpaceConcept<Space>);
    /*!
     * \return the tensorial object built from the components of the view at the
     * given integration point
     * \tparam begin: index of the first component of the tensorial object
     * \tparam ValueType: type of the tensorial object or scalar value returned
     * \param[in] i: offset associated with the integration point
     * \note only available if the components with indices from `begin` to
     * `begin + M - 1` exist, `M` being the number of components of
     * `ValueType`, and for views defined on a linear element space
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(
        const size_type) const  //
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearElementSpaceConcept<Space>);
    /*!
     * \return the tensorial object built from the components of the view at the
     * given integration point
     * \tparam begin: index of the first component of the tensorial object
     * \tparam ValueType: type of the tensorial object or scalar value returned
     * \param[in] e: cell index
     * \param[in] q: quadrature point index
     * \note only available if the components with indices from `begin` to
     * `begin + M - 1` exist, `M` being the number of components of
     * `ValueType`, and for views defined on a linear quadrature space
     */
    template <size_type begin, ScalarOrTensorConcept ValueType>
    [[nodiscard]] constexpr ConstValues<ValueType> get(const size_type,
                                                       const size_type) const
        requires((begin + internals::CompileTimeSize<ValueType>::value <= N) &&
                 LinearQuadratureSpaceConcept<Space>);

  };  // end of StridedCoalescedMemoryAccessCompositeTensorsView

}  // namespace mgis::function

#include "MGIS/Function/TFEL/StridedCoalescedMemoryAccessCompositeTensorsView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_STRIDEDCOALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX \
        */
