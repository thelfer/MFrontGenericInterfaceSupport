/*!
 * \file   MGIS/Function/ModifierConcept.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   05/06/2025
 */

#ifndef LIB_MGIS_FUNCTION_MODIFIERCONCEPT_IXX
#define LIB_MGIS_FUNCTION_MODIFIERCONCEPT_IXX

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType, typename ModifierType>
  constexpr auto operator|(const EvaluatorType& e,
                           ModifierType m) requires(requires {
    { m(e) } -> EvaluatorConcept;
  }) {
    return m(e);
  }  // end of operator |

}// end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_MODIFIERCONCEPT_IXX */
