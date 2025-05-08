/*!
 * \file   MGIS/Function/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_HXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  //  /*!
  //   * \brief assign the evaluator to a function
  //   * \param[in] ctx: execution context
  //   * \param[in] lhs: left hand side
  //   * \param[in] e: right hand side
  //   */
  //   template <size_type N, typename FunctionType, EvaluatorConcept
  //   EvaluatorType>
  //   [[nodiscard]] bool assign(Context&,
  //                             FunctionType&,
  //                             const EvaluatorType&)  //
  //       requires(
  //           (N > 0) &&
  //           ((LinearElementSpaceConcept<std::decay_t<
  //                 decltype(std::declval<EvaluatorType>().getSpace())>>) ||
  //            (LinearQuadratureSpaceConcept<std::decay_t<
  //                 decltype(std::declval<EvaluatorType>().getSpace())>>)) &&
  //           internals::same_decay_type<decltype(std::declval<FunctionType>().getSpace()),
  //                                      decltype(std::declval<EvaluatorType>().getSpace())>);

  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] bool assign(Context&,
                            FunctionType&,
                            EvaluatorType)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(std::declval<EvaluatorType>().getSpace())>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(std::declval<EvaluatorType>().getSpace())>>)) &&
               internals::same_decay_type<
                   decltype(std::declval<FunctionType>().getSpace()),
                   decltype(std::declval<EvaluatorType>().getSpace())>);

}  // end of namespace mgis::function

#include "MGIS/Function/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_HXX */
