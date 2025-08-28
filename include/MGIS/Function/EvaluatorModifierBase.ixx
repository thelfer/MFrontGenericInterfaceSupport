/*!
 * \file   MGIS/Function/EvaluatorModifierBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_IXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_IXX

namespace mgis::function {

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr EvaluatorModifierBase<Child, EvaluatorType>::EvaluatorModifierBase(
      const EvaluatorType& e)
      : evaluator(e) {}  // end of EvaluatorModifierBase

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr EvaluatorModifierBase<Child, EvaluatorType>::EvaluatorModifierBase(
      const EvaluatorModifierBase&) = default;

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr EvaluatorModifierBase<Child, EvaluatorType>::EvaluatorModifierBase(
      EvaluatorModifierBase&&) = default;

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr bool EvaluatorModifierBase<Child, EvaluatorType>::check(
      AbstractErrorHandler& ctx) const {
    return this->evaluator.check(ctx);
  }  // end of check

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr void
  EvaluatorModifierBase<Child, EvaluatorType>::allocateWorkspace() {
    this->evaluator.allocateWorkspace();
  }  // end of allocatWorkspace

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr decltype(auto)
  EvaluatorModifierBase<Child, EvaluatorType>::getSpace() const {
    return internals::disambiguateGetSpace(this->evaluator);
  }  // end of getSpace

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr auto EvaluatorModifierBase<Child, EvaluatorType>::operator()(
      const element_index<Space>& e) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b1) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->evaluator(e));
  }  // end of operator()

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr auto EvaluatorModifierBase<Child, EvaluatorType>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& e) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b2) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->evaluator(wk, e));
  }  // end of operator()

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr auto EvaluatorModifierBase<Child, EvaluatorType>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b3) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->evaluator(e, i));
  }  // end of operator()

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr auto EvaluatorModifierBase<Child, EvaluatorType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b4) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->evaluator(wk, e, i));
  }  // end of operator()

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr decltype(auto) getSpace(
      const EvaluatorModifierBase<Child, EvaluatorType>& e) {
    return e.getSpace();
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_IXX */
