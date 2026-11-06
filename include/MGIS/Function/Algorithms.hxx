/*!
 * \file   MGIS/Function/Algorithms.hxx
 * \brief Declaration of the algorithms which can be applied to the functions
 * and the evaluators of the `MGIS/Function` library, that is `assign` and
 * `scalar_reduce`, and of the helpers used to check if an evaluator can be
 * assigned to a function.
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

#include <version>
// The parallel versions of the algorithms are only built if the library is
// configured with the `enable-parallel-stl-algorithms` option, which defines
// the `MGIS_USE_STL_PARALLEL_ALGORITHMS` macro, and if the compiler provides
// the STL parallel algorithms, as checked using the
// `__cpp_lib_parallel_algorithm` macro defined by the `<version>` header
#ifdef MGIS_USE_STL_PARALLEL_ALGORITHMS
#ifdef __cpp_lib_parallel_algorithm
#define MGIS_HAS_STL_PARALLEL_ALGORITHMS
#endif /* __cpp_lib_parallel_algorithm */
#endif /* MGIS_USE_STL_PARALLEL_ALGORITHMS */

#include <concepts>
#include <optional>
#ifdef MGIS_USE_STL_PARALLEL_ALGORITHMS
#include <execution>
#endif /* MGIS_USE_STL_PARALLEL_ALGORITHMS */
#include "MGIS/Config.hxx"
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"

namespace mgis::function::internals {

  /*!
   * \brief assigns to a span of dynamic size the values of any container which
   * can be iterated over
   * \param[out] lhs: span receiving the values
   * \param[in] rhs: container providing the values
   * \note the number of values of `lhs` must be greater than or equal to the
   * number of values of `rhs`
   */
  constexpr void assign_value(std::span<real>, const auto&);

  /*!
   * \brief assigns to a span of dynamic size the values of a constant span
   */
  constexpr void assign_value(std::span<real>, const std::span<const real>);

  /*!
   * \brief assigns to a span of dynamic size the values of a mutable span
   */
  constexpr void assign_value(std::span<real>, const std::span<real>);

  /*!
   * \brief assigns to a scalar value the first value of a mutable span
   */
  constexpr void assign_value(real&, const std::span<real>&);

  /*!
   * \brief assigns to a scalar value the first value of a constant span
   */
  constexpr void assign_value(real&, const std::span<const real>&);

  /*!
   * \brief assigns to a scalar value the value of a mutable span of one element
   */
  constexpr void assign_value(real&, const std::span<real, 1u>&);

  /*!
   * \brief assigns to a scalar value the value of a constant span of one
   * element
   */
  constexpr void assign_value(real&, const std::span<const real, 1u>&);

  /*!
   * \brief assigns to a scalar value the value of an array of one element
   */
  constexpr void assign_value(real&, const std::array<real, 1u>&);

  /*!
   * \brief assigns a scalar value to the unique element of an array
   */
  constexpr void assign_value(std::array<real, 1u>&, const real);

  /*!
   * \brief assigns the value of an array of one element to an array of one
   * element
   */
  constexpr void assign_value(std::array<real, 1u>&,
                              const std::array<real, 1u>);

  /*!
   * \brief assigns the first value of a mutable span to an array of one element
   */
  constexpr void assign_value(std::array<real, 1u>&, const std::span<real>);

  /*!
   * \brief assigns the first value of a constant span to an array of one
   * element
   */
  constexpr void assign_value(std::array<real, 1u>&,
                              const std::span<const real>);

  /*!
   * \brief assigns the value of a mutable one element span to an array of
   * one element
   */
  constexpr void assign_value(std::array<real, 1u>&, const std::span<real, 1u>);

  /*!
   * \brief assigns the value of a constant span of one element to an array of
   * one element
   */
  constexpr void assign_value(std::array<real, 1u>&,
                              const std::span<const real, 1u>);

  /*!
   * \brief assigns a scalar value to the unique element of a span
   */
  constexpr void assign_value(std::span<real, 1u>, const real);

  /*!
   * \brief assigns a scalar value to the first element of a span of dynamic
   * size
   * \note only the first element of the span is modified
   */
  constexpr void assign_value(std::span<real>, const real);

  /*!
   * \brief assigns the value of a mutable span of one element to the first
   * element of a span of dynamic size
   */
  constexpr void assign_value(std::span<real>, const std::span<real, 1u>);

  /*!
   * \brief assigns the value of a constant span of one element to the first
   * element of a span of dynamic size
   */
  constexpr void assign_value(std::span<real>, const std::span<const real, 1u>);

  /*!
   * \brief assigns the first value of a mutable span to a span of one element
   */
  constexpr void assign_value(std::span<real, 1u>, const std::span<real>);

  /*!
   * \brief assigns the first value of a constant span to a span of one element
   */
  constexpr void assign_value(std::span<real, 1u>, const std::span<const real>);

  /*!
   * \brief assigns the value of a mutable span of one element to a span of one
   * element
   */
  constexpr void assign_value(std::span<real, 1u>, const std::span<real, 1u>);

  /*!
   * \brief assigns the value of a constant span of one element to a span of one
   * element
   */
  constexpr void assign_value(std::span<real, 1u>,
                              const std::span<const real, 1u>);

  /*!
   * \brief assigns the value of an array of one element to a span of one
   * element
   */
  constexpr void assign_value(std::span<real, 1u>, const std::array<real, 1u>&);

  /*!
   * \brief concept stating that an evaluator can be assigned to a function
   *
   * This is the case if both objects are defined on the same space and if the
   * values returned by the evaluator can be assigned to those returned by the
   * function, directly or using one of the `assign_value` helpers.
   *
   * \tparam EvaluatorType: evaluator type
   * \tparam FunctionType: function type
   */
  template <typename EvaluatorType, typename FunctionType>
  concept isEvaluatorAssignableToFunction =
      (EvaluatorConcept<EvaluatorType>)&&  //
      (FunctionConcept<FunctionType>)&&    //
      (std::same_as<function_space<FunctionType>,
                    evaluator_space<EvaluatorType>>)&&  //
      ((requires(function_result<FunctionType> & lhs,
                 const evaluator_result<EvaluatorType>& rhs) { lhs = rhs; }) ||
       (requires(function_result<FunctionType> & lhs,
                 const evaluator_result<EvaluatorType>& rhs) {
         assign_value(lhs, rhs);
       }));
}  // end of namespace mgis::function::internals

namespace mgis::function {

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  /*!
   * \brief concept matching the execution policies of the STL parallel
   * algorithms which can be used by the algorithms of this library
   *
   * \tparam ExecutionPolicy: execution policy considered
   * \note this concept is only defined when the STL parallel algorithms are
   * available
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
   * \brief assigns the values computed by an evaluator to a function using
   * the given execution policy
   *
   * The spaces of the function and the evaluator are compared and the number
   * of components of the values returned by both objects is checked to be the
   * same before the values computed by the evaluator are assigned.
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: type of the function argument, which may be a
   * reference type
   * \tparam EvaluatorType: evaluator type
   * \param[in] ctx: error handler used to report the errors
   * \param[in,out] f: function assigned
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in] e: evaluator providing the values assigned
   * \return true if the assignment succeeded, false otherwise
   * \note the values computed by the evaluator are checked with the `check`
   * function before being assigned
   * \note with the `std::execution::sequenced_policy` policy, the sequential
   * implementation is called, which allows the function to be used in a
   * `constexpr` context
   */
  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            ViewableFunctionArgumentConcept FunctionType,
            EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool assign(AbstractErrorHandler&,
                                      FunctionType&&,
                                      const ExecutionPolicy,
                                      const EvaluatorType)  //
      requires(
          (internals::isEvaluatorAssignableToFunction<
              EvaluatorType,
              std::decay_t<FunctionType>>)&&  //
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>));
  /*!
   * \brief reduces the values computed by a scalar evaluator to a single
   * value using the given execution policy and reduction operator
   *
   * The values computed by the evaluator are checked before being reduced.
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam EvaluatorType: evaluator type
   * \tparam OperatorType: type of the reduction operator, a binary operation
   * taking two scalar values and returning a scalar value
   * \param[in] ctx: error handler used to report the errors
   * \param[in] policy: execution policy used to perform the reduction
   * \param[in] e: evaluator reduced
   * \param[in] op: reduction operator
   * \param[in] initial_value: value the reduction starts from
   * \return the reduced value or an empty `std::optional` if an error is
   * reported
   * \note this function requires the values returned by the evaluator to have
   * a single component
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
   * \brief assigns the values computed by an evaluator to a function, the
   * evaluation being sequential
   *
   * The spaces of the function and the evaluator are compared and the number
   * of components of the values returned by both objects is checked to be the
   * same before the values computed by the evaluator are assigned.
   *
   * \tparam FunctionType: type of the function argument, which may be a
   * reference type
   * \tparam EvaluatorType: evaluator type
   * \param[in] ctx: error handler used to report the errors
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \return true if the assignment succeeded, false otherwise
   * \note the values computed by the evaluator are checked with the `check`
   * function before being assigned. As no execution policy is used, this
   * function can be called in a `constexpr` context
   */
  template <ViewableFunctionArgumentConcept FunctionType,
            EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool assign(AbstractErrorHandler&,
                                      FunctionType&&,
                                      const EvaluatorType)  //
      requires(
          (internals::isEvaluatorAssignableToFunction<
              EvaluatorType,
              std::decay_t<FunctionType>>)&&  //
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>));
  /*!
   * \brief reduces the values computed by a scalar evaluator to a single
   * value using the given reduction operator, the evaluation being sequential
   *
   * The values computed by the evaluator are checked before being reduced.
   *
   * \tparam EvaluatorType: evaluator type
   * \tparam OperatorType: type of the reduction operator, a binary operation
   * taking two scalar values and returning a scalar value
   * \param[in] ctx: error handler used to report the errors
   * \param[in] e: evaluator reduced
   * \param[in] op: reduction operator
   * \param[in] initial_value: value the reduction starts from
   * \return the reduced value or an empty `std::optional` if an error is
   * reported
   * \note this function requires the values returned by the evaluator to have
   * a single component
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
