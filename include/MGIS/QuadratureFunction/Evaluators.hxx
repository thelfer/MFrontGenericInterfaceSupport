/*!
 * \file   MGIS/QuadratureFunction/Evaluators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_HXX

#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/QuadratureFunction/Buffer.hxx"
#include "MGIS/QuadratureFunction/AbstractQuadratureSpace.hxx"
#include "MGIS/QuadratureFunction/QuadratureFunction.hxx"

namespace mgis::quadrature_function {

  template <typename EvaluatorType>
  concept EvaluatorConcept = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> &&
      requires(EvaluatorType& e) {
    e.allocateWorkspace();
  } && requires(const EvaluatorType& e) {
    e.check();
    { e.getQuadratureSpace() } -> std::same_as<const AbstractQuadratureSpace&>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
  } &&((requires(const EvaluatorType& e, size_type i) { e(i); }) ||
       (requires(EvaluatorType & e, size_type i) { e(i); }));

  /*!
   * \brief an evaluator returning the values of an immutable partial
   * quadrature function view as a fixed size span or a scalar \tparam size of
   * the returned value
   */
  template <size_type N>
  struct FixedSizedEvaluator {
    /*!
     * \brief constructor
     */
    FixedSizedEvaluator(const ImmutableQuadratureFunctionView&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const AbstractQuadratureSpace& getQuadratureSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const size_type) const;

   private:
    //! \brief underlying partial quadrature space
    const ImmutableQuadratureFunctionView& function;
    };  // end of FixedSizedEvaluator

    /*!
     * \brief check if the given evaluators have the same partial quadrature
     * space \param[in] e1: first evaluator \param[in] e2: second evaluator
     */
    template <EvaluatorConcept EvaluatorType1, EvaluatorConcept EvaluatorType2>
    void checkMatchingAbstractQuadratureSpaces(const EvaluatorType1&,
                                               const EvaluatorType2&);

  }  // end of namespace mgis::quadrature_function

#include "MGIS/QuadratureFunction/Evaluators.ixx"

#endif /* LIB_MGIS_QUADRATUREFUNCTION_EVALUATORS_HXX */
