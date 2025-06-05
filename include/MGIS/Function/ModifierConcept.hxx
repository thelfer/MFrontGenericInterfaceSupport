/*!
 * \file   ModifierConcept.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   05/06/2025
 */

#ifndef LIB_MGIS_FUNCTION_MODIFIERCONCEPT_HXX
#define LIB_MGIS_FUNCTION_MODIFIERCONCEPT_HXX

#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function{

  /*!
   * \return the evaluator resulting from appling the modifier to the evaluator
   * \param[in] e: evaluator
   * \param[in] m: modifier
   */
  template <EvaluatorConcept EvaluatorType, typename ModifierType>
  constexpr auto operator|(const EvaluatorType& e,
                           ModifierType m) requires(requires {
    { m(e) } -> EvaluatorConcept;
  });

}

#include "MGIS/Function/ModifierConcept.ixx"

#endif /* LIB_MGIS_FUNCTION_MODIFIERCONCEPT_HXX */
