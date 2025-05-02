/*!
 * \file   MGIS/QuadratureFunction/Algorithms.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   23/04/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_ALGORITHMS_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_ALGORITHMS_HXX

#include "MGIS/Config.hxx"
#include "MGIS/QuadratureFunction.hxx"
#include "MGIS/QuadratureFunctionEvaluators.hxx"

namespace mgis::quadrature_function{

  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <size_type N, typename QuadratureFunctionEvaluatorType>
  bool assign(QuadratureFunction&,
              QuadratureFunctionEvaluatorType) requires(N > 0);
  /*!
   * \brief assign the evaluator to a partial quadrature function
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename QuadratureFunctionEvaluatorType>
  bool assign(QuadratureFunction&,
              QuadratureFunctionEvaluatorType);

  /*!
  template <typename ValueType, typename BinaryOperator>
  ValueType reduce(ImmutableQuadratureFunctionView f,
                   const ValueType init,
                   BinaryOperator op) {
    constexpr bool expects_scalar_function =
        requires(const real v, const ValueType& v2) {
      op(v, v2);
    };
    const auto ne =
        f.getQuadratureSpace().getNumberOfIntegrationPoints();
    auto r = init;
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (expects_scalar_function) {
        const auto v = f.getIntegrationPointValue(i);
        r = op(v, r);
      } else {
        const auto v = f.getIntegrationPointValues(i);
        r = op(v, r);
      }
    }
    return r;
  }
   */

}  // end of namespace mgis::quadrature_function

#include "MGIS/QuadratureFunction/Algorithms.ixx"

#endif /* LIB_MGIS_QUADRATUREFUNCTION_ALGORITHMS_HXX */
