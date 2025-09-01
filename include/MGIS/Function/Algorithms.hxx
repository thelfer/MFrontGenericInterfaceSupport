/*!
 * \file   MGIS/Function/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_HXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_HXX

#include <optional>
#include "MGIS/Config.hxx"
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"

namespace mgis::function {

  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool assign(AbstractErrorHandler&,
                                      FunctionType&,
                                      EvaluatorType)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>);

  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] e: evaluator reduced.
   * \param[in] op: reduction operator.
   * \param[in] initial_value: initial value.
   */
  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  [[nodiscard]] std::optional<real> scalar_reduce(
      AbstractErrorHandler&,
      EvaluatorType,
      const OperatorType,
      const real) requires(LinearElementSpaceConcept<std::
                                                         decay_t<decltype(getSpace(
                                                             std::declval<
                                                                 EvaluatorType>()))>>);

}  // end of namespace mgis::function

#include "MGIS/Function/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_HXX */
