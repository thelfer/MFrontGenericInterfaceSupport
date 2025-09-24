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
#include <algorithm>
#include <type_traits>

#ifdef MGIS_HAVE_TFEL
#include "TFEL/Math/General/BasicOperations.hxx"
#include "TFEL/Math/General/UnaryResultType.hxx"
#include "TFEL/Math/General/ResultType.hxx"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/EvaluatorModifierBase.hxx"

namespace mgis::function {

  namespace internals {

    template <typename CallableType, EvaluatorConcept EvaluatorType>
    struct UnaryOperationBase {
      static constexpr size_type getNumberOfComponents() noexcept {
        using UnaryOperationResult =
            std::invoke_result_t<CallableType, evaluator_result<EvaluatorType>>;
        return compile_time_size<UnaryOperationResult>;
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
    constexpr UnaryOperation(const CallableType&, const EvaluatorType&);
    //
    using internals::UnaryOperationBase<CallableType,
                                        EvaluatorType>::getNumberOfComponents;
    //! \brief apply the modifier
    constexpr auto apply(const evaluator_result<EvaluatorType>&) const;

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
    //
    using internals::UnaryOperationBase<CallableType,
                                        EvaluatorType>::getNumberOfComponents;
    //! \brief apply the modifier
    constexpr auto apply(const evaluator_result<EvaluatorType>&) const;
  };

  //! \return the number of components
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation<CallableType, EvaluatorType>&) noexcept;
  //! \return the number of components
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation2<CallableType, EvaluatorType>&) noexcept;

  namespace internals {

    template <typename CallableType>
    struct unary_operation_modifier {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      //
      constexpr unary_operation_modifier(const CallableType&);
      //
      template <typename EvaluatorType>
      constexpr auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));

     private:
      CallableType modifier;
    };

    template <typename CallableType>
    struct unary_operation_modifier2_impl {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      //
      template <typename EvaluatorType>
      constexpr auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));
    };

    template <typename CallableType>
    constexpr auto unary_operation_modifier2(CallableType);

  }  // namespace internals

  template <typename CallableType>
  constexpr auto unary_operation(CallableType&&);

  template <typename CallableType>
  constexpr auto transform(CallableType&&);

}  // end of namespace mgis::function

#include "MGIS/Function/UnaryOperation.ixx"

namespace mgis::function::customization_points {

  template <typename T>
  struct AbsoluteValue;

  inline constexpr auto absolute_value = []<typename T>(const T& v) {
    return AbsoluteValue<T>::exe(v);
  };

  template <>
  struct AbsoluteValue<real> {
    static constexpr real exe(const real& v) noexcept { return v < 0 ? -v : v; }
  };

  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct AbsoluteValue<std::span<const real, N>> {
    static constexpr auto exe(const std::span<const real, N>& v) noexcept {
      auto r = std::array<real, N>{};
      for (std::size_t i = 0; i != N; ++i) {
        r[i] = std::abs(v[i]);
      }
      return r;
    }
  };

  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct AbsoluteValue<std::array<const real, N>> {
    static constexpr auto exe(const std::array<const real, N>& v) noexcept {
      auto r = std::array<real, N>{};
      for (std::size_t i = 0; i != N; ++i) {
        r[i] = std::abs(v[i]);
      }
      return r;
    }
  };

  template <typename T>
  struct MaximumComponent;

  inline constexpr auto maximum_component = []<typename T>(const T& v) {
    return MaximumComponent<T>::exe(v);
  };

  template <>
  struct MaximumComponent<real> {
    static constexpr real exe(const real& v) noexcept { return v; }
  };

  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct MaximumComponent<std::span<const real, N>> {
    static constexpr real exe(const std::span<const real, N>& v) noexcept
        requires((N != 0) && (N != std::dynamic_extent)) {
      if constexpr (N == 1) {
        return v[0];
      } else if constexpr (N == 2) {
        return std::max(v[0], v[1]);
      } else if constexpr (N == 3) {
        return std::max(std::max(v[0], v[1]), v[2]);
      } else if constexpr (N == 4) {
        return std::max(std::max(std::max(v[0], v[1]), v[2]), v[3]);
      } else {
        return *(std::max_element(v.begin(), v.end()));
      }
    }
  };

  template <std::size_t N>
  requires(N > 0) struct MaximumComponent<std::array<real, N>> {
    static constexpr real exe(const std::array<real, N>& v) noexcept {
      return maximum_component(std::span<const real, N>(v.data(), N));
    }
  };

  template <typename T>
  struct MinimumComponent;

  inline constexpr auto minimum_component = []<typename T>(const T& v) {
    return MinimumComponent<T>::exe(v);
  };

  template <>
  struct MinimumComponent<real> {
    static constexpr real exe(const real& v) noexcept { return v; }
  };

  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct MinimumComponent<std::span<const real, N>> {
    static constexpr real exe(const std::span<const real, N>& v) noexcept
        requires((N != 0) && (N != std::dynamic_extent)) {
      if constexpr (N == 1) {
        return v[0];
      } else if constexpr (N == 2) {
        return std::min(v[0], v[1]);
      } else if constexpr (N == 3) {
        return std::min(std::min(v[0], v[1]), v[2]);
      } else if constexpr (N == 4) {
        return std::min(std::min(std::min(v[0], v[1]), v[2]), v[3]);
      } else {
        return *(std::min_element(v.begin(), v.end()));
      }
    }
  };

  template <std::size_t N>
  requires(N > 0) struct MinimumComponent<std::array<real, N>> {
    static constexpr real exe(const std::array<real, N>& v) noexcept {
      return minimum_component(std::span<const real, N>(v.data(), N));
    }
  };

}  // namespace mgis::function::customization_points

namespace mgis::function {

  inline constexpr auto absolute_value = internals::unary_operation_modifier2(
      []<typename ValueType>(const ValueType& v) constexpr {
        return customization_points::absolute_value(v);
      });

  inline constexpr auto maximum_component =
      internals::unary_operation_modifier2(
          []<typename ValueType>(const ValueType& v) constexpr {
            return customization_points::maximum_component(v);
          });

  inline constexpr auto minimum_component =
      internals::unary_operation_modifier2(
          []<typename ValueType>(const ValueType& v) constexpr {
            return customization_points::minimum_component(v);
          });

}  // end of namespace mgis::function

#ifdef MGIS_HAVE_TFEL

namespace mgis::function {

  inline constexpr auto negate = internals::unary_operation_modifier2(
      []<typename OperandType>(const OperandType& a) constexpr
          ->tfel::math::UnaryOperationResult<OperandType,
                                             tfel::math::OpNeg>  //
      requires(compile_time_size<
                   tfel::math::UnaryOperationResult<OperandType,
                                                    tfel::math::OpNeg>> !=
               dynamic_extent) {  //
        return -a;
      });

  constexpr auto multiply_by_scalar(const real s) {
    auto c = [b = s]<typename FirstOperandType>(const FirstOperandType& a)
        -> tfel::math::BinaryOperationResult<FirstOperandType, real,
                                             tfel::math::OpMult>  //
    requires(compile_time_size<tfel::math::BinaryOperationResult<
                 FirstOperandType, real, tfel::math::OpMult>> !=
             dynamic_extent) {
      return a * b;
    };
    return internals::unary_operation_modifier<decltype(c)>(c);
  }

  constexpr auto divide_by_scalar(const real s) {
    auto c = [b = s]<typename FirstOperandType>(const FirstOperandType& a)
        -> tfel::math::BinaryOperationResult<FirstOperandType, real,
                                             tfel::math::OpDiv>  //
    requires(compile_time_size<tfel::math::BinaryOperationResult<
                 FirstOperandType, real, tfel::math::OpDiv>> !=
             dynamic_extent) {
      return a / b;
    };
    return internals::unary_operation_modifier<decltype(c)>(c);
  }  // end of divide_by_scalar

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#endif /* LIB_MGIS_FUNCTION_UNARYOPERATION_HXX */
