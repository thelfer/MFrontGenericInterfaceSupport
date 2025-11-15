/*!
 * \file   MGIS/Function/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/04/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_HXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_HXX

#include <concepts>
#include <optional>
#include <execution>
#include "MGIS/Config.hxx"
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"

#ifdef MGIS_USE_STL_PARALLEL_ALGORITHMS
#ifdef __cpp_lib_parallel_algorithm
#define MGIS_HAS_STL_PARALLEL_ALGORITHMS
#endif /* __cpp_lib_parallel_algorithm */
#endif /* MGIS_USE_STL_PARALLEL_ALGORITHMS */

namespace mgis::function {

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  /*!
   * \brief concept matching one of the supported execution policies
   */
  template <typename ExecutionPolicy>
  concept ExecutionPolicyConceptConcept =
      std::same_as<ExecutionPolicy, std::execution::sequenced_policy> ||
      std::same_as<ExecutionPolicy, std::execution::unsequenced_policy> ||
      std::same_as<ExecutionPolicy, std::execution::parallel_policy> ||
      std::same_as<ExecutionPolicy,
                   std::execution::parallel_unsequenced_policy>;

#endif

#ifndef _MSC_VER

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS
  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] policy: execution policy
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool assign(AbstractErrorHandler&,
                                      FunctionType&,
                                      const ExecutionPolicy,
                                      const EvaluatorType)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>);
  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] policy: execution policy
   * \param[in] e: evaluator reduced.
   * \param[in] op: reduction operator.
   * \param[in] initial_value: initial value.
   */
  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  [[nodiscard]] constexpr std::optional<real> scalar_reduce(
      AbstractErrorHandler&,
      const ExecutionPolicy,
      const EvaluatorType,
      const OperatorType,
      const real) requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>);
#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool assign(AbstractErrorHandler&,
                                      FunctionType&,
                                      const EvaluatorType)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>);
  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] e: evaluator reduced.
   * \param[in] op: reduction operator.
   * \param[in] initial_value: initial value.
   */
  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  [[nodiscard]] constexpr std::optional<real> scalar_reduce(
      AbstractErrorHandler&,
      const EvaluatorType,
      const OperatorType,
      const real) requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>);

#endif /* MSC_VER */

}  // end of namespace mgis::function

#include "MGIS/Function/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_HXX */
