/*!
 * \file   MGIS/Function/Evaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATOR_HXX
#define LIB_MGIS_FUNCTION_EVALUATOR_HXX

#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief check if the given evaluators shares the same space
   *
   * \param[in] ctx: context
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   */
  bool checkMatchingSpaces(Context&,
                           const EvaluatorConcept auto&,
                           const EvaluatorConcept auto&);

  /*!
   * \return the evaluator resulting from appling the modifier to the evaluator
   * \param[in] e: evaluator
   * \param[in] m: modifier
   */
  template <EvaluatorConcept EvaluatorType, typename ModifierType>
  auto operator|(EvaluatorType,
                 ModifierType)  //
      requires(requires(EvaluatorType e1, ModifierType m1) {
        { m1(e1) } -> EvaluatorConcept;
      });

}  // end of namespace mgis::function

#include "MGIS/Function/Evaluator.ixx"
#include "MGIS/Function/UnaryOperation.hxx"
#include "MGIS/Function/BinaryOperation.hxx"

#endif /* LIB_MGIS_FUNCTION_EVALUATOR_HXX */
