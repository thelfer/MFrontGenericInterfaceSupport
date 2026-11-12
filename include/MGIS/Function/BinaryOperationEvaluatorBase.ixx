/*!
 * \file   MGIS/Function/BinaryOperationEvaluatorBase.ixx
 * \brief  Implementation of the `BinaryOperationEvaluatorBase` class
 * declared in `MGIS/Function/BinaryOperationEvaluatorBase.hxx`.
 * \author Thomas Helfer
 * \date   07/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX

namespace mgis::function {

  //! \brief constructor, storing a copy of each of the two evaluators
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

  //! \brief default copy constructor
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

  //! \brief default move constructor
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

  /*!
   * \brief implementation of the `check` member function
   *
   * The spaces of the two evaluators are compared first, then the values
   * computed by each evaluator are checked. The first failing check stops the
   * process and registers the associated error message in the error handler.
   */
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
    if (!internals::disambiguateCheck(ctx, this->first_evaluator)) {
      return false;
    }
    if (!internals::disambiguateCheck(ctx, this->second_evaluator)) {
      return false;
    }
    return true;
  }  // end of check

  // the space of the first evaluator is returned: the requirements state
  // that the two evaluators are defined on the same space and the `check`
  // function verifies that their spaces are equivalent
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

  // the call operators downcast the object to the child class and pass the
  // values computed by the two evaluators to its `apply` member function
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

  //! \brief forwards the call to the `getSpace` member function
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr decltype(auto) getSpace(
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>& e) {
    return e.getSpace();
  }  // end of getSpace

  //! \brief forwards the call to the `check` member function
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr bool check(
      AbstractErrorHandler& eh,
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>& e) {
    return e.check(eh);
  }  // end of check

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_IXX */
