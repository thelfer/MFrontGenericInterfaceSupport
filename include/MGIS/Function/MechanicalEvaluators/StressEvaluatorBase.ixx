/*!
 * \file   StressEvaluatorBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX
#define LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"

namespace mgis::function {

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::
      StressEvaluatorBase(const StressEvaluatorType& e)
      : stress_evaluator(e) {}  // end of StressEvaluatorBase

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::
      StressEvaluatorBase(const StressEvaluatorBase&) = default;

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::
      StressEvaluatorBase(StressEvaluatorBase&&) = default;

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  void StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::check()
      const {
    using namespace tfel::math;
    this->stress_evaluator.check();
    const auto nc = this->stress_evaluator.getNumberOfComponents();
    if constexpr (isSymmetric) {
      raise_if(nc != StensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    } else {
      raise_if(nc != TensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    }
  }  // end of check

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  void StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::
      allocateWorkspace() {
    this->stress_evaluator.allocatWorkspace();
  }  // end of allocatWorkspace

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  const auto&
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::getSpace()
      const {
    return this->stress_evaluator.getSpace();
  }  // end of getSpace

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  auto
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::operator()(
      const element_index<Space>& e) const
      requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->stress_evaluator(e));
  }  // end of operator()

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  auto
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& e) const
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->stress_evaluator(wk, e));
  }  // end of operator()

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  auto
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::operator()(
      const cell_index<Space> e, const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>)) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->stress_evaluator(e, i));
  }  // end of operator()

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  auto
  StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space> e,
      const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    const auto& child = static_cast<const Child&>(*this);
    return child.apply(this->stress_evaluator(wk, e, i));
  }  // end of operator()

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX */
