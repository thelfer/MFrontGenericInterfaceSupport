
/*!
 * \file   MGIS/Function/TFEL/TensorialFunction.ixx
 * \brief  File which contains the implementation of the `TensorialFunction`
 * class declared in `MGIS/Function/TFEL/TensorialFunction.hxx`.
 * \author Thomas Helfer
 * \date   28/08/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_IXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_IXX

namespace mgis::function {

  // the constructor checking the preconditions delegates to the one allowing
  // the preconditions to be skipped
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorialFunction<Space, TensorType>::TensorialFunction(
          const Space& s)
      : TensorialFunction(preconditions_check, s) {
  }  // end of TensorialFunction

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool TensorialFunction<Space, TensorType>::checkPreconditions(
          AbstractErrorHandler&, const Space&) {
    return true;
  }

  // the function owns the values, stored by the `Function` base class, and the
  // view, which shares them, is built from this base class, so it refers to the
  // values of the function it belongs to
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr TensorialFunction<Space, TensorType>::TensorialFunction(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s)
      : PreconditionsChecker<TensorialFunction>(pcheck, s),
        Function<Space,
                 tfel::math::getUnderlyingArrayMinimalSize<
                     typename TensorType::indexing_policy>()>(s),
        TensorView<Function<Space,
                            tfel::math::getUnderlyingArrayMinimalSize<
                                typename TensorType::indexing_policy>()>,
                   TensorType>(
            pcheck,
            static_cast<Function<Space,
                                 tfel::math::getUnderlyingArrayMinimalSize<
                                     typename TensorType::indexing_policy>()>&>(
                *this)) {}  // end of TensorialFunction

  // the view is rebound to the values copied from the source function; the
  // preconditions are not checked as the view is built from the values of the
  // function itself
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorialFunction<Space, TensorType>::TensorialFunction(
          const TensorialFunction& src)
      : PreconditionsChecker<TensorialFunction>(src),
        UnderlyingFunctionType(static_cast<const UnderlyingFunctionType&>(src)),
        UnderlyingTensorView(no_precondition_check,
                             static_cast<UnderlyingFunctionType&>(*this)) {
  }  // end of TensorialFunction

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorialFunction<Space, TensorType>::TensorialFunction(
          TensorialFunction&& src)
      : PreconditionsChecker<TensorialFunction>(src),
        UnderlyingFunctionType(static_cast<UnderlyingFunctionType&&>(src)),
        UnderlyingTensorView(no_precondition_check,
                             static_cast<UnderlyingFunctionType&>(*this)) {
  }  // end of TensorialFunction

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorView<
          Function<Space,
                   tfel::math::getUnderlyingArrayMinimalSize<
                       typename TensorType::indexing_policy>()>,
          TensorType>&  //
      TensorialFunction<Space, TensorType>::view() & MGIS_NOEXCEPT {
    return *this;
  }  // end of view

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr const
      TensorView<Function<Space,
                          tfel::math::getUnderlyingArrayMinimalSize<
                              typename TensorType::indexing_policy>()>,
                 TensorType>&  //
      TensorialFunction<Space, TensorType>::view() const& MGIS_NOEXCEPT {
    return *this;
  }  // end of view

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorialFunction<Space, TensorType>::~TensorialFunction() =
          default;

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr decltype(auto) getSpace(
      const TensorialFunction<Space, TensorType>& f) {
    return f.getSpace();
  }

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(TensorialFunction<Space, TensorType>& f) {
    return f.view();
  }  // end of view

  // a read-only `FunctionView` of the values owned by the function is created
  // and wrapped in a `TensorModifier` so that its values are returned as
  // read-only tensors viewing these values
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(const TensorialFunction<Space, TensorType>& f) {
    constexpr auto N = tfel::math::getUnderlyingArrayMinimalSize<
        typename TensorType::indexing_policy>();
    using FunctionViewType =
        FunctionView<Space, simple_data_layout_description<N>, false>;
    return TensorModifier<TensorType, FunctionViewType>(
        FunctionViewType(f.getSpace(), f.data()));
  }  // end of view

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorialFunction<Space, TensorType>& f) MGIS_NOEXCEPT {
    return f.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_IXX */
