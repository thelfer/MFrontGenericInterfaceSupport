/*!
 * \file   MGIS/Function/TensorialFunction.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   28/08/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORIALFUNCTION_IXX
#define LIB_MGIS_FUNCTION_TENSORIALFUNCTION_IXX

#ifdef MGIS_HAVE_TFEL

namespace mgis::function {

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

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorView<
          Function<Space,
                   tfel::math::getUnderlyingArrayMinimalSize<
                       typename TensorType::indexing_policy>()>,
          TensorType>&  //
      TensorialFunction<Space, TensorType>::view() noexcept {
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
      TensorialFunction<Space, TensorType>::view() const noexcept {
    return *this;
  }  // end of view

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr TensorialFunction<Space, TensorType>::~TensorialFunction() =
          default;

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr decltype(auto)
      getSpace(const TensorialFunction<Space, TensorType>& f) {
    return f.getSpace();
  }

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(const TensorialFunction<Space, TensorType>& f){
    return f.view();
  } // end of view

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorialFunction<Space, TensorType>& f) noexcept {
    return f.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#endif /* LIB_MGIS_FUNCTION_TENSORIALFUNCTION_HXX */
