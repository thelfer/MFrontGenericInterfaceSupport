/*!
 * \file   MGIS/Function/EvaluatorConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   21/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX
#define LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX

namespace mgis::function {

  constexpr bool checkMatchingSpaces(AbstractErrorHandler& ctx,
                                     const EvaluatorConcept auto& e1,
                                     const EvaluatorConcept auto& e2) {
    const auto& qspace1 = getSpace(e1);
    const auto& qspace2 = getSpace(e2);
    if (!areEquivalent(qspace1, qspace2)) {
      return ctx.registerErrorMessage("unmatched quadrature spaces");
    }
    return true;
  }  // end of checkMatchingSpaces

}  // end of namespace mgis::function

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
  constexpr mgis::size_type disambiguateGetNumberOfComponents(
      const EvaluatorType& e) {
    return getNumberOfComponents(e);
  }  // end of disambiguateGetNumberOfComponents

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX */
