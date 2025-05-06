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
  concept EvaluatorConceptBase = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> &&
      requires(EvaluatorType& e) {
    e.allocateWorkspace();
  } && requires(const EvaluatorType& e) {
    e.check();
    FunctionalSpaceConcept<decltype(e.getSpace())>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
  };

  template <typename EvaluatorType>
  concept LinearEvaluatorConcept = EvaluatorConceptBase<EvaluatorType> &&
      LinearSpaceConcept<decltype(std::declval<EvaluatorType>().getSpace())> &&
      ((requires(const EvaluatorType& e, size_type i) { e(i); }) ||
       (requires(EvaluatorType & e, size_type i) { e(i); }));

  template <typename EvaluatorType>
  concept QuadratureEvaluatorConcept = EvaluatorConceptBase<EvaluatorType> &&
      QuadratureSpaceConcept<
          decltype(std::declval<EvaluatorType>().getSpace())> &&
      ((requires(const EvaluatorType& e, size_type n, size_type i) {
         e(n, i);
       }) ||
       (requires(EvaluatorType & e, size_type n, size_type i) { e(n, i); }));

  template <typename EvaluatorType>
  concept EvaluatorConcept = LinearEvaluatorConcept<std::decay_t<EvaluatorType>> ||
      QuadratureEvaluatorConcept<std::decay_t<EvaluatorType>>;

  template <typename EvaluatorType>
  concept LinearQuadratureEvaluatorConcept =
      LinearEvaluatorConcept<std::decay_t<EvaluatorType>> &&
      QuadratureEvaluatorConcept<std::decay_t<EvaluatorType>>;

  /*!
   * \brief an evaluator returning the values of an immutable partial
   * quadrature function view as a fixed size span or a scalar
   *
   * \tparam Space: functional space
   * \tparam N: size of the returned value
   */
  template <FunctionalSpaceConcept Space, size_type N>
  struct FixedSizeEvaluator {
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] values: function
     */
    static bool checkPreconditions(
        const ImmutableFunctionView<Space, {}>&) noexcept;
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    FixedSizeEvaluator(const ImmutableFunctionView<Space, {}>&);
    //! \brief perform consistency checks
    bool check() const noexcept;
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
    auto operator()(const typename SpaceTraits<Space>::element_index_type) const
        requires(LinearSpaceConcept<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(
        const typename SpaceTraits<Space>::CellWorkspace&,
        const typename SpaceTraits<Space>::cell_index_type,
        const typename SpaceTraits<Space>::quadrature_point_index_type) const
        requires(LinearQuadratureSpaceConcept<Space>);

   private:
    //! \brief underlying partial quadrature space
    const ImmutableFunctionView<Space, {}>& function;
  };  // end of FixedSizeEvaluator

  /*!
   * \brief check if the given evaluators have the same partial quadrature
   * space
   *
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   */
  void checkMatchingAbstractSpaces(const EvaluatorConceptBase auto&,
                                   const EvaluatorConceptBase auto&);

}  // end of namespace mgis::function

#include "MGIS/Function/Evaluators.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORS_HXX */
