/*!
 * \file   MGIS/Function/FunctionConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX
#define LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX

#include "MGIS/Context.hxx"
#include "MGIS/Function/Algorithms.hxx"

namespace mgis::function::internals {

  template <FunctionConcept FunctionType>
  decltype(auto) disambiguateGetSpace(const FunctionType& f)  //
      requires(!EvaluatorConcept<FunctionType>) {
    return getSpace(f);
  }

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType, FunctionConcept FunctionType>
  bool operator|(EvaluatorType e, FunctionType& f) requires(
      internals::same_decay_type<
          decltype(getSpace(std::declval<EvaluatorType>())),
          decltype(getSpace(std::declval<FunctionType>()))>) {
    Context ctx;
    return assign(ctx, f, e);
  }  // end of operator |

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX */
