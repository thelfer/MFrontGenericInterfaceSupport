/*!
 * \file   MGIS/Function/EvaluatorModifierConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   05/06/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType,
            EvaluatorModifierConcept EvaluatorModifierType>
  constexpr auto operator|(const EvaluatorType& e, EvaluatorModifierType m) {
    return m(e);
  }  // end of operator |

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX */
