/*!
 * \file   MGIS/Function/EvaluatorConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   21/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX
#define LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX

namespace mgis::function::internals {

  template <EvaluatorConcept EvaluatorType>
  constexpr decltype(auto) disambiguateGetSpace(const EvaluatorType& e) {
    return getSpace(e);
  }

  template <EvaluatorConcept EvaluatorType>
  constexpr bool disambiguateCheck(AbstractErrorHandler& ctx,
                                   const EvaluatorType& e) {
    return check(ctx, e);
  }

  template <EvaluatorConcept EvaluatorType>
  constexpr void disambiguateAllocateWorkspace(EvaluatorType& e) {
    allocateWorkspace(e);
  }  // end of disambiguateAllocateWorkspace

  template <EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type disambiguateGetNumberOfComponents(
      const EvaluatorType& e) {
    return getNumberOfComponents(e);
  }  // end of disambiguateGetNumberOfComponents

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX */
