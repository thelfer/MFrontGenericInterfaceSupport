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

  template <size_type N>
  FixedSizedEvaluator<N>::FixedSizedEvaluator(
      const ImmutableFunctionView& values)
      : function(values) {}

  template <size_type N>
  void FixedSizedEvaluator<N>::check() const {
    raise_if(this->function.getNumberOfComponents() != N,
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }

  template <size_type N>
  void FixedSizedEvaluator<N>::allocateWorkspace() {}

  template <size_type N>
  const AbstractSpace& FixedSizedEvaluator<N>::getSpace() const {
    return this->function.getSpace();
  }

  template <size_type N>
  constexpr size_type FixedSizedEvaluator<N>::getNumberOfComponents()
      const noexcept {
    return N;
  }

  template <size_type N>
  auto FixedSizedEvaluator<N>::operator()(const size_type i) const {
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
