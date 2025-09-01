/*!
 * \file   MGIS/Function/BinaryOperationEvaluatorBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX

namespace mgis::function {

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr BinaryOperationEvaluatorBase<Child,
                                             FirstEvaluatorType,
                                             SecondEvaluatorType>::
          BinaryOperationEvaluatorBase(const FirstEvaluatorType& e1,
                                       const SecondEvaluatorType& e2)
      : first_evaluator(e1),
        second_evaluator(e2) {}  // end of BinaryOperationEvaluatorBase

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr BinaryOperationEvaluatorBase<Child,
                                             FirstEvaluatorType,
                                             SecondEvaluatorType>::
          BinaryOperationEvaluatorBase(const BinaryOperationEvaluatorBase&) =
              default;

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr BinaryOperationEvaluatorBase<Child,
                                             FirstEvaluatorType,
                                             SecondEvaluatorType>::
          BinaryOperationEvaluatorBase(BinaryOperationEvaluatorBase&&) =
              default;

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr bool BinaryOperationEvaluatorBase<
          Child,
          FirstEvaluatorType,
          SecondEvaluatorType>::check(AbstractErrorHandler& ctx) const {
    if (!checkMatchingSpaces(ctx, this->first_evaluator,
                             this->second_evaluator)) {
      return false;
    }
    if (!this->first_evaluator.check(ctx)) {
      return false;
    }
    if (!this->second_evaluator.check(ctx)) {
      return false;
    }
    return true;
  }  // end of check

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr void BinaryOperationEvaluatorBase<
          Child,
          FirstEvaluatorType,
          SecondEvaluatorType>::allocateWorkspace() {
    internals::disambiguateAllocateWorkspace(this->first_evaluator);
    internals::disambiguateAllocateWorkspace(this->second_evaluator);
  }  // end of allocatWorkspace

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr decltype(auto)
          BinaryOperationEvaluatorBase<Child,
                                       FirstEvaluatorType,
                                       SecondEvaluatorType>::getSpace() const {
    return internals::disambiguateGetSpace(this->first_evaluator);
  }  // end of getSpace

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr auto BinaryOperationEvaluatorBase<
          Child,
          FirstEvaluatorType,
          SecondEvaluatorType>::operator()(const element_index<Space>& e) const
      requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b1) &&
               (internals::EvaluatorResultQuery<SecondEvaluatorType>::b1)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(e), this->second_evaluator(e));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr auto BinaryOperationEvaluatorBase<
          Child,
          FirstEvaluatorType,
          SecondEvaluatorType>::operator()(const element_workspace<Space>& wk,
                                           const element_index<Space>& e) const
      requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b2) &&
               (internals::EvaluatorResultQuery<SecondEvaluatorType>::b2)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(wk, e),
                       this->second_evaluator(wk, e));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr auto BinaryOperationEvaluatorBase<Child,
                                                  FirstEvaluatorType,
                                                  SecondEvaluatorType>::
      operator()(const cell_index<Space> e,
                 const quadrature_point_index<Space> i) const
      requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b3) &&
               (internals::EvaluatorResultQuery<SecondEvaluatorType>::b3)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(e, i),
                       this->second_evaluator(e, i));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      constexpr auto BinaryOperationEvaluatorBase<Child,
                                                  FirstEvaluatorType,
                                                  SecondEvaluatorType>::
      operator()(const cell_workspace<Space>& wk,
                 const cell_index<Space> e,
                 const quadrature_point_index<Space> i) const
      requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b4) &&
               (internals::EvaluatorResultQuery<SecondEvaluatorType>::b4)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(wk, e, i),
                       this->second_evaluator(wk, e, i));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr decltype(auto) getSpace(
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>& e) {
    return e.getSpace();
  }  // end of getSpace

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr void allocateWorkspace(
      BinaryOperationEvaluatorBase<Child,
                                   FirstEvaluatorType,
                                   SecondEvaluatorType>& e) {
    return e.allocateWorkspace();
  }  // end of allocateWorkspace

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX */
