/*!
 * \file   MGIS/Function/Evaluators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORS_HXX
#define LIB_MGIS_FUNCTION_EVALUATORS_HXX

#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Function/Buffer.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  template <typename EvaluatorType>
  concept EvaluatorConcept = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> &&
      requires(EvaluatorType& e) {
    e.allocateWorkspace();
  } && requires(const EvaluatorType& e) {
    e.check();
    FunctionalSpaceConcept<decltype(e.getSpace())>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
  } &&(LinearSpaceConcept<decltype(std::declval<EvaluatorType>().getSpace())>
           ? ((requires(const EvaluatorType& e, size_type i) { e(i); }) ||
              (requires(EvaluatorType & e, size_type i) { e(i); }))
           : true) &&
      (QuadratureSpaceConcept<
           decltype(std::declval<EvaluatorType>().getSpace())>
           ? ((requires(const EvaluatorType& e, size_type n, size_type i) {
                e(n, i);
              }) ||
              (requires(EvaluatorType & e, size_type n, size_type i) {
                e(n, i);
              }))
           : true);

  /*!
   * \brief an evaluator returning the values of an immutable partial
   * quadrature function view as a fixed size span or a scalar \tparam size of
   * the returned value
   */
  template <FunctionalSpaceConcept Space, size_type N>
  struct FixedSizedEvaluator {
    /*!
     * \brief constructor
     */
    FixedSizedEvaluator(const ImmutableFunctionView<Space, {}>&) noexcept;
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const Space& getSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const size_type) const;

   private:
    //! \brief underlying partial quadrature space
    const ImmutableFunctionView<Space, {}>& function;
  };  // end of FixedSizedEvaluator

  /*!
   * \brief check if the given evaluators have the same partial quadrature
   * space \param[in] e1: first evaluator \param[in] e2: second evaluator
   */
  template <EvaluatorConcept EvaluatorType1, EvaluatorConcept EvaluatorType2>
  void checkMatchingAbstractSpaces(const EvaluatorType1&,
                                   const EvaluatorType2&);

}  // end of namespace mgis::function

#include "MGIS/Function/Evaluators.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORS_HXX */
