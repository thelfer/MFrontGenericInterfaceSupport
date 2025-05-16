/*!
 * \file   MGIS/Function/BinaryOperation.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATION_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATION_HXX

#ifdef MGIS_HAVE_TFEL
#include "TFEL/Math/General/BasicOperations.hxx"
#include "TFEL/Math/General/ResultType.hxx"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/UnaryOperation.hxx"
#include "MGIS/Function/BinaryOperationEvaluatorBase.hxx"

namespace mgis::function {

  namespace internals {

    template <typename CallableType,
              EvaluatorConcept FirstEvaluatorType,
              EvaluatorConcept SecondEvaluatorType>
    struct BinaryOperationModifierBase {
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
  inline constexpr auto BinaryOperationModifierRequirements =
      (std::is_copy_constructible_v<CallableType>)&&(
          std::invocable<CallableType,
                         evaluator_result<FirstEvaluatorType>,
                         evaluator_result<SecondEvaluatorType>>);

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifierRequirements<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType>)  //
      struct BinaryOperationModifier
      : BinaryOperationEvaluatorBase<
            BinaryOperationModifier<CallableType,
                                    FirstEvaluatorType,
                                    SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationModifierBase<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType> {
    /*!
     * \brief constructor
     * \param[in] c: callable type
     * \param[in] e1: evaluator associated with the first argument
     * \param[in] e2: evaluator associated with the second argument
     */
    BinaryOperationModifier(const CallableType&,
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
  inline constexpr auto BinaryOperationModifier2Requirements =
      (std::is_trivially_default_constructible_v<CallableType>)&&(
          std::invocable<CallableType,
                         evaluator_result<FirstEvaluatorType>,
                         evaluator_result<SecondEvaluatorType>>);

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifier2Requirements<CallableType,
                                                FirstEvaluatorType,
                                                SecondEvaluatorType>)  //
      struct BinaryOperationModifier2
      : BinaryOperationEvaluatorBase<
            BinaryOperationModifier2<CallableType,
                                     FirstEvaluatorType,
                                     SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationModifierBase<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType> {
    // inheriting constructor
    using BinaryOperationEvaluatorBase<
        BinaryOperationModifier2,
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
      binary_operation_modifier(const CallableType&);
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
    struct binary_operation_modifier2_impl {
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

    template <typename CallableType>
    constexpr auto binary_operation_modifier2(CallableType);

  }  // namespace internals

  template <typename CallableType, typename SecondEvaluatorType>
  auto binary_operation(CallableType&&,
                        SecondEvaluatorType&&)  //
      requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);

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

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperation.ixx"

#ifdef MGIS_HAVE_TFEL

namespace mgis::function {

  inline constexpr auto add = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<FirstOperandType,
                                               SecondOperandType,
                                               tfel::math::OpPlus>  //
      requires(compile_time_size<
                   tfel::math::BinaryOperationResult<FirstOperandType,
                                                     SecondOperandType,
                                                     tfel::math::OpPlus>> !=
               dynamic_extent) { return a + b; });

  inline constexpr auto substract = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<FirstOperandType,
                                               SecondOperandType,
                                               tfel::math::OpMinus>  //
      requires(compile_time_size<
                   tfel::math::BinaryOperationResult<FirstOperandType,
                                                     SecondOperandType,
                                                     tfel::math::OpMinus>> !=
               dynamic_extent) { return a - b; });

  inline constexpr auto mean_value = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<
              tfel::math::BinaryOperationResult<FirstOperandType,
                                                SecondOperandType,
                                                tfel::math::OpMinus>,
              real,
              tfel::math::OpDiv>  //
      requires(compile_time_size<tfel::math::BinaryOperationResult<
                   tfel::math::BinaryOperationResult<FirstOperandType,
                                                     SecondOperandType,
                                                     tfel::math::OpMinus>,
                   real,
                   tfel::math::OpDiv>> != dynamic_extent) {
        return (a + b) / 2;
      });

  inline constexpr auto multiply = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<FirstOperandType,
                                               SecondOperandType,
                                               tfel::math::OpMult>  //
      requires(compile_time_size<
                   tfel::math::BinaryOperationResult<FirstOperandType,
                                                     SecondOperandType,
                                                     tfel::math::OpMult>> !=
               dynamic_extent) { return a * b; });

  inline constexpr auto divide = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<FirstOperandType,
                                               SecondOperandType,
                                               tfel::math::OpDiv>  //
      requires(compile_time_size<
                   tfel::math::BinaryOperationResult<FirstOperandType,
                                                     SecondOperandType,
                                                     tfel::math::OpDiv>> !=
               dynamic_extent) { return a / b; });

  inline constexpr auto inner_product = internals::binary_operation_modifier2(
      []<typename FirstOperandType, typename SecondOperandType>(
          const FirstOperandType& a, const SecondOperandType& b)
          -> tfel::math::BinaryOperationResult<FirstOperandType,
                                               SecondOperandType,
                                               tfel::math::OpDotProduct>  //
      requires(
          compile_time_size<
              tfel::math::BinaryOperationResult<FirstOperandType,
                                                SecondOperandType,
                                                tfel::math::OpDotProduct>> !=
          dynamic_extent) { return a | b; });

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#endif /* LIB_MGIS_FUNCTIONBINARYOPERATION_HXX */
