/*!
 * \file   MGIS/Function/TransformEvaluatorModifier.hxx
 * \brief  This file declares the TransformEvaluatorModifier class and the
 * transform function
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_HXX

#include <span>
#include "MGIS/Function/EvaluatorModifierBase.hxx"

namespace mgis::function {

  namespace internals {

    template <EvaluatorConcept EvaluatorType, typename CallableType>
    struct TransformEvaluatorModifierBase {
      static constexpr size_type getNumberOfComponents() noexcept {
        using result_type =
            std::invoke_result_t<CallableType, evaluator_result<EvaluatorType>>;
        return compile_time_size<result_type>;
      }
      static_assert(getNumberOfComponents() != dynamic_extent);
    };

  }  // end of namespace internals

  /*!
   * \brief a base class for evaluators modifying a stress tensor
   * \tparam EvaluatorType: evaluator of the stress
   * \tparam CallableType: type of the modifier
   */
  template <EvaluatorConcept EvaluatorType, typename CallableType>
  requires(std::invocable<CallableType,
                          evaluator_result<EvaluatorType>>)  //
      struct TransformEvaluatorModifier
      : EvaluatorModifierBase<
            TransformEvaluatorModifier<EvaluatorType, CallableType>,
            EvaluatorType>,
        internals::TransformEvaluatorModifierBase<EvaluatorType, CallableType> {
    /*!
     * \brief constructor
     * \param[in] e: modified evaluator
     * \param[in] c: callable type
     */
    TransformEvaluatorModifier(const EvaluatorType&, const CallableType&);
    //! \brief apply the modifier
    auto apply(const evaluator_result<EvaluatorType>&) const;

   private:
    CallableType modifier;
  };

  namespace internals {

    template <typename CallableType>
    struct transform_modifier {
      //
      transform_modifier(CallableType&&);
      //
      template <typename EvaluatorType>
      auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType, evaluator_result<EvaluatorType>>));

     private:
      CallableType modifier;
    };

  }  // namespace internals

  template <typename CallableType>
  auto transform(CallableType&&);

}  // end of namespace mgis::function

#include "MGIS/Function/TransformEvaluatorModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_HXX */
