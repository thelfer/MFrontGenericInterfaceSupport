/*!
 * \file   MGIS/Function/Evaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATOR_HXX
#define LIB_MGIS_FUNCTION_EVALUATOR_HXX

#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Function/Space.hxx"

namespace mgis::function {

  template <typename EvaluatorType>
  concept EvaluatorConcept = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> &&
      requires(EvaluatorType& e) {
    e.allocateWorkspace();
  } && 
    (ElementSpaceConcept<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>&&
        hasElementWorkspace<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>
        ? requires(const EvaluatorType& e) {
     {e2(std::declval<element_workspace<std::decay_t<decltype(e.getSpace())>>>(),
         std::declval<element_index<std::decay_t<decltype(e.getSpace())>>>())};
    }
:true) &&
     (ElementSpaceConcept<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>> &&
         (!hasElementWorkspace<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>)
         ? requires(const EvaluatorType& e) {
      {e(std::declval<element_index<std::decay_t<decltype(e.getSpace())>>>())};
 }
 :true)&&
  (QuadratureSpaceConcept<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>&&
         hasCellWorkspace<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>
         ? requires(const EvaluatorType& e) {
      {e(std::declval<cell_workspace<std::decay_t<decltype(e.getSpace())>>>(),
          std::declval<cell_index<std::decay_t<decltype(e.getSpace())>>>(),
          std::declval<quadrature_point_index<std::decay_t<decltype(e.getSpace())>>>())};
 }
 :true)&&
  (QuadratureSpaceConcept<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>> &&
         (!hasCellWorkspace<std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>>)
         ? requires(const EvaluatorType& e) {
        {e(std::declval<cell_index<std::decay_t<decltype(e.getSpace())>>>(),
          std::declval<quadrature_point_index<std::decay_t<decltype(e.getSpace())>>>())};
 }
 :true)&&
requires(const EvaluatorType& e, Context& ctx) {
    FunctionalSpaceConcept<std::decay_t<decltype(e.getSpace())>>;
    { e.check(ctx) } -> std::same_as<bool>;
    { e.getNumberOfComponents() } -> std::same_as<size_type>;
  };

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

#include "MGIS/Function/Evaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATOR_HXX */
