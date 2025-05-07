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
#include "MGIS/Context.hxx"
#include "MGIS/Function/Buffer.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  template <typename EvaluatorType>
  concept EvaluatorConcept = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> &&
      requires(EvaluatorType& e) {
    e.allocateWorkspace();
  } && requires(const EvaluatorType& e, Context& ctx) {
    FunctionalSpaceConcept<std::decay_t<decltype(e.getSpace())>>;
    { e.check(ctx) } -> std::same_as<bool>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
    (ElementSpaceConcept<std::decay_t<decltype(e.getSpace())>>&&
        hasElementWorkspace<std::decay_t<decltype(e.getSpace())>>
        ? requires(const EvaluatorType& e2) {
      e2(std::declval<element_workspace<std::decay_t<decltype(e.getSpace())>>>(),
         std::declval<element_index<std::decay_t<decltype(e.getSpace())>>>());
    }
:true);
    (ElementSpaceConcept<std::decay_t<decltype(e.getSpace())>> &&
        (!hasElementWorkspace<std::decay_t<decltype(e.getSpace())>>)
        ? requires(const EvaluatorType& e2) {
      e2(std::declval<element_index<std::decay_t<decltype(e.getSpace())>>>());
}
:true);
    (QuadratureSpaceConcept<std::decay_t<decltype(e.getSpace())>>&&
        hasCellWorkspace<std::decay_t<decltype(e.getSpace())>>
        ? requires(const EvaluatorType& e2) {
      e2(std::declval<cell_workspace<std::decay_t<decltype(e.getSpace())>>>(),
         std::declval<cell_index<std::decay_t<decltype(e.getSpace())>>>(),
         std::declval<quadrature_point_index<std::decay_t<decltype(e.getSpace())>>>());
}
:true);
    (QuadratureSpaceConcept<std::decay_t<decltype(e.getSpace())>> &&
        (!hasCellWorkspace<std::decay_t<decltype(e.getSpace())>>)
        ? requires(const EvaluatorType& e2) {
      e2(std::declval<cell_index<std::decay_t<decltype(e.getSpace())>>>(),
         std::declval<quadrature_point_index<std::decay_t<decltype(e.getSpace())>>>());
}
:true);
  };

  /*!
   * \brief an evaluator returning the values of an immutable
   * function view as a fixed size span or a scalar
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
    bool check(Context&) const noexcept;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying  space
    const Space& getSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_index<Space>&) const
        requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_workspace<Space>&,
                    const element_index<Space>&) const
        requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_workspace<Space>&,
                    const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>);

   private:
    //! \brief underlying function
    const ImmutableFunctionView<Space, {}>& function;
  };  // end of FixedSizeEvaluator

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

}  // end of namespace mgis::function

#include "MGIS/Function/Evaluators.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORS_HXX */
