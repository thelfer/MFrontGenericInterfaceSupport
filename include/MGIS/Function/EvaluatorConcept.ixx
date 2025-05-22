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
  decltype(auto) disambiguateGetSpace(const EvaluatorType& e) {
    return getSpace(e);
  }

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_EVALUATORCONCEPT_IXX */
