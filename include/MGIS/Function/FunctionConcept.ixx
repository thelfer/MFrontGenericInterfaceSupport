/*!
 * \file   MGIS/Function/FunctionConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX
#define LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX

#include "MGIS/Function/Algorithms.hxx"

namespace mgis::function {

  /*!
   * \brief assign an evaluator to a mutable function view
   * \param[in] e: evaluator
   * \param[in] f: function
   */
  template <EvaluatorConcept EvaluatorType, FunctionConcept FunctionType>
  bool operator|(EvaluatorType e, FunctionType& f) requires(
      internals::same_decay_type<
          decltype(std::declval<EvaluatorType>().getSpace()),
          decltype(std::declval<FunctionType>().getSpace())>) {
    Context ctx;
    return assign(ctx, f, e);
  }  // end of operator |

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX */
