/*!
 * \file   MGIS/Function/Evaluators.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORS_IXX
#define LIB_MGIS_FUNCTION_EVALUATORS_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N>
  bool FixedSizeEvaluator<Space, N>::checkPreconditions(
      const ImmutableFunctionView<Space, {}>& values) noexcept {
    return values.getNumberOfComponents() == N;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N>
  FixedSizeEvaluator<Space, N>::FixedSizeEvaluator(
      const ImmutableFunctionView<Space, {}>& values)
      : function(values) {
    raise_if(!checkPreconditions(values),
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }  // end of FixedSizeEvaluator

  template <FunctionalSpaceConcept Space, size_type N>
  bool FixedSizeEvaluator<Space, N>::check() const noexcept {
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  void FixedSizeEvaluator<Space, N>::allocateWorkspace() {}

  template <FunctionalSpaceConcept Space, size_type N>
  const Space& FixedSizeEvaluator<Space, N>::getSpace() const {
    return this->function.getSpace();
  }

  template <FunctionalSpaceConcept Space, size_type N>
  constexpr size_type FixedSizeEvaluator<Space, N>::getNumberOfComponents()
      const noexcept {
    return N;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  auto FixedSizeEvaluator<Space, N>::operator()(
      const typename SpaceTraits<Space>::element_index_type i) const
      requires(LinearSpaceConcept<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N>
  auto FixedSizeEvaluator<Space, N>::operator()(
      const typename SpaceTraits<Space>::CellWorkspace&,
      const typename SpaceTraits<Space>::cell_index_type e,
      const typename SpaceTraits<Space>::quadrature_point_index_type i) const
      requires(LinearQuadratureSpaceConcept<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(e, i);
    } else {
      return this->function.template getValues<N>(e, i);
    }
  }

  void checkMatchingAbstractSpaces(const EvaluatorConceptBase auto& e1,
                                   const EvaluatorConceptBase auto& e2) {
    const auto& qspace1 = e1.getSpace();
    const auto& qspace2 = e2.getSpace();
    raise_if(&qspace1 != &qspace2, "unmatched quadrature spaces");
  }  // end of checkMatchingAbstractSpaces

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORS_IXX */
