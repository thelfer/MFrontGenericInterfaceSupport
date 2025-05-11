/*!
 * \file   MGIS/Function/BinaryOperationEvaluator.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   11/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATOR_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/TransformEvaluatorModifier.hxx"
#include "MGIS/Function/BinaryOperationEvaluatorBase.hxx"

namespace mgis::function {

  namespace internals {

    template <typename CallableType,
              EvaluatorConcept FirstEvaluatorType,
              EvaluatorConcept SecondEvaluatorType>
    struct BinaryOperationEvaluatorModifierBase {
      static constexpr size_type getNumberOfComponents() noexcept {
        using result_type =
            std::invoke_result_t<CallableType,
                                 evaluator_result<FirstEvaluatorType>,
                                 evaluator_result<SecondEvaluatorType>>;
        return compile_time_size<result_type>;
      }
      static_assert(getNumberOfComponents() != dynamic_extent);
    };

  }  // end of namespace internals

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<FirstEvaluatorType>,
                     evaluator_result<SecondEvaluatorType>>))  //
      struct BinaryOperationEvaluatorModifier
      : BinaryOperationEvaluatorBase<
            BinaryOperationEvaluatorModifier<CallableType,
                                             FirstEvaluatorType,
                                             SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationEvaluatorModifierBase<CallableType,
                                                        FirstEvaluatorType,
                                                        SecondEvaluatorType> {
    /*!
     * \brief constructor
     * \param[in] c: callable type
     * \param[in] e1: evaluator associated with the first argument
     * \param[in] e2: evaluator associated with the second argument
     */
    BinaryOperationEvaluatorModifier(const CallableType&,
                                     const FirstEvaluatorType&,
                                     const SecondEvaluatorType&);
    //! \brief apply the modifier
    auto apply(const evaluator_result<FirstEvaluatorType>&,
               const evaluator_result<SecondEvaluatorType>&) const;

   private:
    CallableType modifier;
  };

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires((std::is_trivially_default_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<FirstEvaluatorType>,
                     evaluator_result<SecondEvaluatorType>>))  //
      struct BinaryOperationEvaluatorModifier2
      : BinaryOperationEvaluatorBase<
            BinaryOperationEvaluatorModifier2<CallableType,
                                              FirstEvaluatorType,
                                              SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationEvaluatorModifierBase<CallableType,
                                                        FirstEvaluatorType,
                                                        SecondEvaluatorType> {
    // inheriting constructor
    using BinaryOperationEvaluatorBase<
        BinaryOperationEvaluatorModifier2,
        FirstEvaluatorType,
        SecondEvaluatorType>::BinaryOperationEvaluatorBase;
    //! \brief apply the modifier
    auto apply(const evaluator_result<FirstEvaluatorType>&,
               const evaluator_result<SecondEvaluatorType>&) const;
  };

  namespace internals {

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_copy_constructible_v<
             CallableType>) struct BinaryOperatorCurrying {
      //
      BinaryOperatorCurrying(const CallableType&, const SecondEvaluatorType&);
      //
      template <typename FirstEvaluatorType>
      auto operator()(FirstEvaluatorType&&) const requires(
          (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
          (std::invocable<CallableType,
                          evaluator_result<std::decay_t<FirstEvaluatorType>>,
                          evaluator_result<SecondEvaluatorType>>));

     private:
      CallableType modifier;
      SecondEvaluatorType e2;
    };

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_trivially_default_constructible_v<
             CallableType>) struct BinaryOperatorCurrying2 {
      //
      BinaryOperatorCurrying2(const SecondEvaluatorType&);
      //
      template <typename FirstEvaluatorType>
      auto operator()(FirstEvaluatorType&&) const requires(
          (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
          (std::invocable<CallableType,
                          evaluator_result<std::decay_t<FirstEvaluatorType>>,
                          evaluator_result<SecondEvaluatorType>>));

     private:
      SecondEvaluatorType e2;
    };

    template <typename CallableType>
    struct binary_operation_modifier {
      //
      binary_operation_modifier(CallableType&&);
      //
      template <typename SecondEvaluatorType>
      auto operator()(SecondEvaluatorType&&) const
          requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);
      //
      template <typename FirstEvaluatorType, typename SecondEvaluatorType>
      auto operator()(FirstEvaluatorType&&, SecondEvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                   (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                   (std::invocable<
                       CallableType,
                       evaluator_result<std::decay_t<FirstEvaluatorType>>,
                       evaluator_result<std::decay_t<SecondEvaluatorType>>>));

     private:
      CallableType modifier;
    };

    template <typename CallableType>
    struct binary_operation_modifier2 {
      //
      template <typename SecondEvaluatorType>
      auto operator()(SecondEvaluatorType&&) const
          requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);
      //
      template <typename FirstEvaluatorType, typename SecondEvaluatorType>
      auto operator()(FirstEvaluatorType&&, SecondEvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                   (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                   (std::invocable<
                       CallableType,
                       evaluator_result<std::decay_t<FirstEvaluatorType>>,
                       evaluator_result<std::decay_t<SecondEvaluatorType>>>));
    };

  }  // namespace internals

  template <typename CallableType,
            typename FirstEvaluatorType,
            typename SecondEvaluatorType>
  auto binary_operation(CallableType&&,
                        FirstEvaluatorType&&,
                        SecondEvaluatorType&&)                          //
      requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
               (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
               (std::invocable<
                   std::decay_t<CallableType>,
                   evaluator_result<std::decay_t<FirstEvaluatorType>>,
                   evaluator_result<std::decay_t<SecondEvaluatorType>>>));

  //   template <typename CallableType, typename SecondEvaluatorType>
  //   auto binary_operation(CallableType&&, SecondEvaluatorType&&) requires(
  //       EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperationEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTIONBINARYOPERATIONEVALUATOR_HXX */
