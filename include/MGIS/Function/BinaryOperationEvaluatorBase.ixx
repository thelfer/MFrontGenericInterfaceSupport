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
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
      BinaryOperationEvaluatorBase(const FirstEvaluatorType& e1,
                                   const SecondEvaluatorType& e2)
      : first_evaluator(e1),
        second_evaluator(e2) {}  // end of BinaryOperationEvaluatorBase

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
      BinaryOperationEvaluatorBase(const BinaryOperationEvaluatorBase&) =
          default;

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
      BinaryOperationEvaluatorBase(BinaryOperationEvaluatorBase&&) = default;

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  bool BinaryOperationEvaluatorBase<Child,
                                    FirstEvaluatorType,
                                    SecondEvaluatorType>::check(Context& ctx)
      const {
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
  void BinaryOperationEvaluatorBase<Child,
                                    FirstEvaluatorType,
                                    SecondEvaluatorType>::allocateWorkspace() {
    this->first_evaluator.allocatWorkspace();
    this->second_evaluator.allocatWorkspace();
  }  // end of allocatWorkspace

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  const auto& BinaryOperationEvaluatorBase<Child,
                                           FirstEvaluatorType,
                                           SecondEvaluatorType>::getSpace()
      const {
    return this->first_evaluator.getSpace();
  }  // end of getSpace

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
  operator()(const element_index<Space>& e) const
      requires(BinaryOperationEvaluatorBase::isElementEvaluator &&
               !(hasElementWorkspace<Space>)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(e), this->second_evaluator(e));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
  operator()(const element_workspace<Space>& wk,
             const element_index<Space>& e) const
      requires(BinaryOperationEvaluatorBase::isElementEvaluator&&
                   hasElementWorkspace<Space>) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(wk, e),
                       this->second_evaluator(wk, e));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
  operator()(const cell_index<Space> e,
             const quadrature_point_index<Space> i) const
      requires(BinaryOperationEvaluatorBase::isQuadratureEvaluator &&
               (!hasCellWorkspace<Space>)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(e, i),
                       this->second_evaluator(e, i));
  }  // end of operator()

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto
  BinaryOperationEvaluatorBase<Child, FirstEvaluatorType, SecondEvaluatorType>::
  operator()(const cell_workspace<Space>& wk,
             const cell_index<Space> e,
             const quadrature_point_index<Space> i) const
      requires(BinaryOperationEvaluatorBase::isQuadratureEvaluator&&
                   hasCellWorkspace<Space>) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->first_evaluator(wk, e, i),
                       this->second_evaluator(wk, e, i));
  }  // end of operator()

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX */
