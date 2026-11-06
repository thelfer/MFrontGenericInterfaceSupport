/*!
 * \file   MGIS/Function/TFEL/CoalescedMemoryAccessCompositeTensorsView.hxx
 * \brief  File which defines the `CoalescedMemoryAccessCompositeTensorsView`
 * class, a view handling several tensorial objects of possibly different types
 * sharing the same set of components, whose values are stored in separate
 * contiguous memory regions.
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

#ifndef LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX

#include <tuple>
#include <concepts>
#include <type_traits>
#include "TFEL/Math/Array/CoalescedView.hxx"
#include "MGIS/Function/CoalescedMemoryAccessFunctionViewBase.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief traits defining the type used by the
   * `CoalescedMemoryAccessCompositeTensorsView` class to return a mutable
   * value of the given type
   *
   * The primary template returns a coalesced view of the tensorial type, whose
   * components are stored in separate contiguous memory regions. It is
   * specialized below for scalar values, which are returned by reference.
   *
   * \tparam ValueType: type of the value to be returned, which must match the
   * `ScalarOrTensorConcept` concept
   */
  template <ScalarOrTensorConcept ValueType>
  struct CoalescedMemoryAccessCompositeTensorsViewMutableValue {
    //! \brief type of the mutable value returned
    using type = tfel::math::CoalescedView<ValueType>;
  };

  /*!
   * \brief specialization for scalar values, which are returned as mutable
   * references to the value of the corresponding component
   */
  template <>
  struct CoalescedMemoryAccessCompositeTensorsViewMutableValue<real> {
    //! \brief type of the mutable value returned
    using type = real&;
  };

  /*!
   * \brief traits defining the type used by the
   * `CoalescedMemoryAccessCompositeTensorsView` class to return a read-only
   * value of the given type
   *
   * The primary template returns a coalesced view of the const tensorial type,
   * whose components are stored in separate contiguous memory regions. It is
   * specialized below for scalar values, which are returned by const reference.
   *
   * \tparam ValueType: type of the value to be returned, which must match the
   * `ScalarOrTensorConcept` concept
   */
  template <ScalarOrTensorConcept ValueType>
  struct CoalescedMemoryAccessCompositeTensorsViewConstValue {
    //! \brief type of the read-only value returned
    using type = tfel::math::CoalescedView<const ValueType>;
  };

  /*!
   * \brief specialization for scalar values, which are returned as read-only
   * references to the value of the corresponding component
   */
  template <>
  struct CoalescedMemoryAccessCompositeTensorsViewConstValue<real> {
    //! \brief type of the read-only value returned
    using type = const real&;
  };

  /*!
   * \brief a coalesced memory access view handling several tensorial objects
   * of possibly different types sharing the same set of components
   *
   * The view is built either from the views of the `N` scalar components or
   * from a space and an array of values split into such views, as described by
   * the `CoalescedMemoryAccessFunctionViewBase` class. The values of each
   * component are stored contiguously for the different points of the space, so
   * the values of the different components at a given integration point are not
   * contiguous in memory.
   *
   * The `get` member function returns, for a given integration point, the
   * tensorial object built from the components with indices from `begin` to
   * `begin + M - 1`, `M` being the number of components of the queried type.
   * Several tensorial objects and scalar values may thus be extracted from the
   * same view, provided that their ranges of components do not exceed the
   * `N` available components, as stated by the `requires`-clause.
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
      struct CoalescedMemoryAccessCompositeTensorsView
      : CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable> {
    //
    /*!
     * \brief type returned by the mutable `get` member functions
     *
     * \tparam ValueType: type of the tensorial object or scalar value returned
     */
    template <ScalarOrTensorConcept ValueType>
    using MutableValues =
        typename CoalescedMemoryAccessCompositeTensorsViewMutableValue<
            ValueType>::type;
    //
    /*!
     * \brief type returned by the constant `get` member functions
     *
     * \tparam ValueType: type of the tensorial object or scalar value returned
     */
    template <ScalarOrTensorConcept ValueType>
    using ConstValues =
        typename CoalescedMemoryAccessCompositeTensorsViewConstValue<
            ValueType>::type;
    /*!
     * \brief constructors inherited from the
     * `CoalescedMemoryAccessFunctionViewBase` base class
     */
    using CoalescedMemoryAccessFunctionViewBase<Space, N, is_mutable>::
        CoalescedMemoryAccessFunctionViewBase;
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

  };  // end of CoalescedMemoryAccessCompositeTensorsView

}  // namespace mgis::function

#include "MGIS/Function/TFEL/CoalescedMemoryAccessCompositeTensorsView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_COALESCEDMEMORYACCESSCOMPOSITETENSORVIEW_HXX \
        */
