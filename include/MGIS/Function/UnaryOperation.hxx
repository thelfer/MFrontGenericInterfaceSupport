/*!
 * \file   MGIS/Function/UnaryOperation.hxx
 * \brief  This file declares the UnaryOperation class and the
 * transform function
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_UNARYOPERATION_HXX
#define LIB_MGIS_FUNCTION_UNARYOPERATION_HXX

#include <span>
#include <type_traits>

#ifdef MGIS_HAVE_TFEL
#include "TFEL/Math/General/BasicOperations.hxx"
#include "TFEL/Math/General/UnaryResultType.hxx"
#include "TFEL/Math/General/ResultType.hxx"
#endif MGIS_HAVE_TFEL

#include "MGIS/Function/EvaluatorModifierBase.hxx"

namespace mgis::function {

  namespace internals {

    template <typename CallableType, EvaluatorConcept EvaluatorType>
    struct UnaryOperationBase {
      static constexpr size_type getNumberOfComponents() noexcept {
        using BinaryOperationResult =
            std::invoke_result_t<CallableType, evaluator_result<EvaluatorType>>;
        return compile_time_size<BinaryOperationResult>;
      }
      static_assert(getNumberOfComponents() != dynamic_extent);
    };

  }  // end of namespace internals

  /*!
   * \brief a base class for evaluators modifying a stress tensor
   * \tparam CallableType: type of the modifier
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      struct UnaryOperation
      : EvaluatorModifierBase<UnaryOperation<CallableType, EvaluatorType>,
                              EvaluatorType>,
        internals::UnaryOperationBase<CallableType, EvaluatorType> {
    /*!
     * \brief constructor
     * \param[in] c: callable type
     * \param[in] e: modified evaluator
     */
    UnaryOperation(const CallableType&, const EvaluatorType&);
    //! \brief apply the modifier
    auto apply(const evaluator_result<EvaluatorType>&) const;

   private:
    CallableType modifier;
  };

  /*!
   * \brief a base class for evaluators modifying a stress tensor
   * \tparam CallableType: type of the modifier
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_trivially_default_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      struct UnaryOperation2
      : EvaluatorModifierBase<UnaryOperation2<CallableType, EvaluatorType>,
                              EvaluatorType>,
        internals::UnaryOperationBase<CallableType, EvaluatorType> {
    // inheriting constructor
    using EvaluatorModifierBase<UnaryOperation2,
                                EvaluatorType>::EvaluatorModifierBase;
    //! \brief apply the modifier
    auto apply(const evaluator_result<EvaluatorType>&) const;
  };

  namespace internals {

    template <typename CallableType>
    struct unary_operation_modifier {
      //
      unary_operation_modifier(const CallableType&);
      //
      template <typename EvaluatorType>
      auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));

     private:
      CallableType modifier;
    };

    template <typename CallableType>
    struct unary_operation_modifier2_impl {
      template <typename EvaluatorType>
      auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));
    };

    template <typename CallableType>
    constexpr auto unary_operation_modifier2(CallableType);

  }  // namespace internals

  template <typename CallableType>
  auto unary_operation(CallableType&&);

  template <typename CallableType>
  auto transform(CallableType&&);

}  // end of namespace mgis::function

#include "MGIS/Function/UnaryOperation.ixx"

#ifdef MGIS_HAVE_TFEL

namespace mgis::function {

  inline constexpr auto negate = internals::unary_operation_modifier2(
      []<typename OperandType>(const OperandType& a)
          -> tfel::math::UnaryOperationResult<OperandType,
                                                     tfel::math::OpNeg>  //
      requires(
          compile_time_size<
              tfel::math::UnaryOperationResult<OperandType,
                                                      tfel::math::OpNeg>> !=
          dynamic_extent) {  //
        return -a;
      });

  inline auto multiply_by_scalar(const real s) {
    auto c = [b = s]<typename FirstOperandType>(const FirstOperandType& a)
        -> tfel::math::BinaryOperationResult<FirstOperandType, real,
                                             tfel::math::OpMult>  //
    requires(compile_time_size<tfel::math::BinaryOperationResult<
                 FirstOperandType, real, tfel::math::OpMult>> !=
             dynamic_extent) {
      return a * b;
    };
    return internals::unary_operation_modifier(c);
  }

  inline auto divide_by_scalar(const real s) {
    auto c = [b = s]<typename FirstOperandType>(const FirstOperandType& a)
        -> tfel::math::BinaryOperationResult<FirstOperandType, real,
                                             tfel::math::OpDiv>  //
    requires(compile_time_size<tfel::math::BinaryOperationResult<
                 FirstOperandType, real, tfel::math::OpDiv>> !=
             dynamic_extent) {
      return a / b;
    };
    return internals::unary_operation_modifier(c);
  }  // end of divide_by_scalar

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#endif /* LIB_MGIS_FUNCTION_UNARYOPERATION_HXX */
