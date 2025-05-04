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
  FixedSizedEvaluator<Space, N>::FixedSizedEvaluator(
      const ImmutableFunctionView<Space, {}>& values) noexcept
      : function(values) {}  // end of FixedSizedEvaluator

  template <FunctionalSpaceConcept Space, size_type N>
  void FixedSizedEvaluator<Space, N>::check() const {
    raise_if(this->function.getNumberOfComponents() != N,
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }

  template <FunctionalSpaceConcept Space, size_type N>
  void FixedSizedEvaluator<Space, N>::allocateWorkspace() {}

  template <FunctionalSpaceConcept Space, size_type N>
  const Space& FixedSizedEvaluator<Space, N>::getSpace() const {
    return this->function.getSpace();
  }

  template <FunctionalSpaceConcept Space, size_type N>
  constexpr size_type FixedSizedEvaluator<Space, N>::getNumberOfComponents()
      const noexcept {
    return N;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  auto FixedSizedEvaluator<Space, N>::operator()(const size_type i) const {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <EvaluatorConcept EvaluatorType1, EvaluatorConcept EvaluatorType2>
  void checkMatchingAbstractSpaces(const EvaluatorType1& e1,
                                   const EvaluatorType2& e2) {
    const auto& qspace1 = e1.getSpace();
    const auto& qspace2 = e2.getSpace();
    raise_if(&qspace1 != &qspace2, "unmatched quadrature spaces");
  }  // end of checkMatchingAbstractSpaces

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORS_IXX */
