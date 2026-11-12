/*!
 * \file   MGIS/Function/UnaryOperation.ixx
 * \brief  Implementation of the `UnaryOperation` and `UnaryOperation2` classes,
 * of the modifiers and of the helpers declared in
 * `MGIS/Function/UnaryOperation.hxx`.
 * \author Thomas Helfer
 * \date   09/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_UNARYOPERATION_IXX
#define LIB_MGIS_FUNCTION_UNARYOPERATION_IXX

namespace mgis::function {

  // the callable object is copied in the object, the underlying
  // evaluator is handled by the `EvaluatorModifierBase` base class
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      constexpr UnaryOperation<CallableType, EvaluatorType>::UnaryOperation(
          const CallableType& c, const EvaluatorType& e)
      : EvaluatorModifierBase<UnaryOperation<CallableType, EvaluatorType>,
                              EvaluatorType>(e),
        modifier(c) {}  // end of UnaryOperation

  //! \brief calls the stored callable object with the given values
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      MGIS_HOST_DEVICE
      constexpr auto UnaryOperation<CallableType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    return this->modifier(values);
  }  // end of apply

  // the callable object is not stored in the object: a default
  // constructed object is used each time the values are computed, which is
  // allowed since `CallableType` is trivially default constructible
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_trivially_default_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      MGIS_HOST_DEVICE
      constexpr auto UnaryOperation2<CallableType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    auto c = CallableType{};
    return c(values);
  }  // end of apply

  //! \brief forwards the call to the `getNumberOfComponents` member
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation<CallableType, EvaluatorType>& e) MGIS_NOEXCEPT {
    return e.getNumberOfComponents();
  }

  //! \brief forwards the call to the `getNumberOfComponents` member
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation2<CallableType, EvaluatorType>& e) MGIS_NOEXCEPT {
    return e.getNumberOfComponents();
  }

  namespace internals {

    //! \brief constructor, storing a copy of the callable object
    template <typename CallableType>
    constexpr unary_operation_modifier<CallableType>::unary_operation_modifier(
        const CallableType& c)
        : modifier(c) {}

    //! \brief builds a `UnaryOperation` from the given evaluator
    template <typename CallableType>
    template <typename EvaluatorType>
    constexpr auto unary_operation_modifier<CallableType>::operator()(
        EvaluatorType&& e) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return UnaryOperation<CallableType, std::decay_t<EvaluatorType>>(
          this->modifier, std::forward<EvaluatorType>(e));
    }  // end of operator()

    //! \brief builds a `UnaryOperation2` from the given evaluator
    template <typename CallableType>
    template <typename EvaluatorType>
    MGIS_HOST_DEVICE constexpr auto
    unary_operation_modifier2_impl<CallableType>::operator()(EvaluatorType&& e)
        const requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return UnaryOperation2<CallableType, std::decay_t<EvaluatorType>>(
          std::forward<EvaluatorType>(e));
    }  // end of operator()

    //! \brief the argument is only used to deduce the type of the
    //!        unary operation
    template <typename CallableType>
    constexpr auto unary_operation_modifier2(CallableType) {
      return unary_operation_modifier2_impl<CallableType>{};
    }  // end of unary_operation_modifier2

  }  // namespace internals

  //! \brief creates a modifier storing a copy of the given callable
  template <typename CallableType>
  constexpr auto unary_operation(CallableType&& c) {
    return internals::unary_operation_modifier<std::decay_t<CallableType>>(
        std::forward<CallableType>(c));
  }  // end of unary_operation

  //! \brief applies the modifier to the given evaluator using `operator|`
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr auto unary_operation(CallableType&& c, const EvaluatorType& e) {
    return e | unary_operation(c);
  }  // end of unary_operation

  //! \brief alias of the single argument `unary_operation` function
  template <typename CallableType>
  constexpr auto transform(CallableType&& c) {
    return internals::unary_operation_modifier<std::decay_t<CallableType>>(
        std::forward<CallableType>(c));
  }  // end of transform

  //! \brief alias of the two arguments `unary_operation` function
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr auto transform(CallableType&& c, const EvaluatorType& e) {
    return e | transform(c);
  }  // end of transform

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_UNARYOPERATION_IXX */
