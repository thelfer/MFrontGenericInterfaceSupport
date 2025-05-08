/*!
 * \file   MGIS/Function/Evaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATOR_IXX
#define LIB_MGIS_FUNCTION_EVALUATOR_IXX

namespace mgis::function {

  bool checkMatchingSpaces(Context& ctx,
                           const EvaluatorConcept auto& e1,
                           const EvaluatorConcept auto& e2) {
    const auto& qspace1 = e1.getSpace();
    const auto& qspace2 = e2.getSpace();
    if (&qspace1 != &qspace2) {
      return ctx.registerErrorMessage("unmatched quadrature spaces");
    }
    return true;
  }  // end of checkMatchingSpaces

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATOR_IXX */
