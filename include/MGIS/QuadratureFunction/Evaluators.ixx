/*!
 * \file   MGIS/QuadratureFunction/Evaluators.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_IXX
#define LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_IXX

#include "MGIS/Raise.hxx"

namespace mgis::quadrature_function {

  template <size_type N>
  FixedSizedEvaluator<N>::FixedSizedEvaluator(
      const ImmutableQuadratureFunctionView& values)
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
  const AbstractQuadratureSpace& FixedSizedEvaluator<N>::getQuadratureSpace()
      const {
    return this->function.getQuadratureSpace();
  }

  template <size_type N>
  constexpr size_type FixedSizedEvaluator<N>::getNumberOfComponents()
      const noexcept {
    return N;
  }

  template <size_type N>
  auto FixedSizedEvaluator<N>::operator()(const size_type i) const {
    if constexpr (N == 1) {
      return this->function.getIntegrationPointValue(i);
    } else {
      return this->function.template getIntegrationPointValues<N>(i);
    }
  }

  template <EvaluatorConcept EvaluatorType1, EvaluatorConcept EvaluatorType2>
  void checkMatchingAbstractQuadratureSpaces(const EvaluatorType1& e1,
                                             const EvaluatorType2& e2) {
    const auto& qspace1 = e1.getQuadratureSpace();
    const auto& qspace2 = e2.getQuadratureSpace();
    raise_if(&qspace1 != &qspace2, "unmatched quadrature spaces");
  }  // end of checkMatchingAbstractQuadratureSpaces

}  // end of namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_IXX */
