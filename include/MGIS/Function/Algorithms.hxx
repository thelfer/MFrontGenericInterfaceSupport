/*!
 * \file   MGIS/Function/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_HXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function.hxx"
#include "MGIS/FunctionEvaluators.hxx"

namespace mgis::function {

  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <size_type N, typename FunctionEvaluatorType>
  bool assign(Function&, FunctionEvaluatorType) requires(N > 0);
  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename FunctionEvaluatorType>
  bool assign(Function&, FunctionEvaluatorType);

  /*!
  template <typename ValueType, typename BinaryOperator>
  ValueType reduce(ImmutableFunctionView f,
                   const ValueType init,
                   BinaryOperator op) {
    constexpr bool expects_scalar_function =
        requires(const real v, const ValueType& v2) {
      op(v, v2);
    };
    const auto ne =
        f.getSpace().getSpaceSize();
    auto r = init;
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (expects_scalar_function) {
        const auto v = f.getValue(i);
        r = op(v, r);
      } else {
        const auto v = f.getValues(i);
        r = op(v, r);
      }
    }
    return r;
  }
   */

}  // end of namespace mgis::function

#include "MGIS/Function/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_HXX */
