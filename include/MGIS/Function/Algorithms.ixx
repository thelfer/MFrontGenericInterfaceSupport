/*!
 * \file   MGIS/Function/Algorithm.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_IXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_IXX

#include <ranges>
#include <numeric>
#include <iterator>
#include <algorithm>

namespace mgis::function::algorithm {

  /*!
   * \brief copy N values
   * \tparam N: number of values to be copied
   * \tparam InputIterator: input iterator type
   * \tparam OutputIterator: output iterator type
   * \param[in] p: iterator to the beginning of the values to be copied
   * \param[in] pe: iterator past the end of the values to be copied
   * \param[in] po: iterator to the output values
   */
  template <size_type N, typename InputIterator, typename OutputIterator>
  constexpr void copy(const InputIterator p,
                      const InputIterator pe,
                      OutputIterator po) requires(N > 0) {
    if constexpr ((std::random_access_iterator<InputIterator>)&&  //
                  (std::random_access_iterator<OutputIterator>)) {
      if constexpr (N > 9) {
        std::copy(p, pe, po);
      } else if constexpr (N == 1) {
        *po = *p;
      } else if constexpr (N == 2) {
        po[0] = p[0];
        po[1] = p[1];
      } else if constexpr (N == 3) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
      } else if constexpr (N == 4) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
        po[3] = p[3];
      } else {
        copy<N - 1>(++p, pe, ++po);
      }
    } else {
      std::copy(p, pe, po);
    }
  }  // end of copy

}  // end of namespace mgis::function::algorithm

namespace mgis::function::internals {

  constexpr void assign_value(auto& lhs, const auto& rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  constexpr void assign_value(std::span<real> lhs, const auto& rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  constexpr void assign_value(std::span<real> lhs,
                              const std::span<const real> rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  constexpr void assign_value(std::span<real> lhs, const std::span<real> rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  constexpr void assign_value(real& lhs, const std::span<real>& rhs) {
    lhs = rhs[0];
  }

  constexpr void assign_value(real& lhs, const std::span<const real>& rhs) {
    lhs = rhs[0];
  }

  constexpr void assign_value(real& lhs, const std::span<real, 1u>& rhs) {
    lhs = rhs[0];
  }

  constexpr void assign_value(real& lhs, const std::span<const real, 1u>& rhs) {
    lhs = rhs[0];
  }

  constexpr void assign_value(real& lhs, const std::array<real, 1u>& rhs) {
    lhs = rhs[0];
  }

  constexpr void assign_value(std::array<real, 1u>& lhs, const real rhs) {
    lhs[0] = rhs;
  }

  constexpr void assign_value(std::array<real, 1u>& lhs,
                              const std::array<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::array<real, 1u>& lhs,
                              const std::span<real> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::array<real, 1u>& lhs,
                              const std::span<const real> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::array<real, 1u>& lhs,
                              const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::array<real, 1u>& lhs,
                              const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs, const real rhs) {
    lhs[0] = rhs;
  }

  constexpr void assign_value(std::span<real> lhs, const real rhs) {
    lhs[0] = rhs;
  }

  constexpr void assign_value(std::span<real> lhs,
                              const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real> lhs,
                              const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs,
                              const std::span<real> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs,
                              const std::span<const real> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs,
                              const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs,
                              const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  constexpr void assign_value(std::span<real, 1u> lhs,
                              const std::array<real, 1u>& rhs) {
    lhs[0] = rhs[0];
  }

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_scalar_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      std::for_each(iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) { f(i) = e(i); });
    } else {
      std::for_each(policy, iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) { f(i) = e(i); });
    }
  }  // end of assign_scalar_impl

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_scalar_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      std::for_each(iranges.begin(), iranges.end(),
                    [&space, &f, e](const space_size_type i) {
                      const auto& wk = getElementWorkspace(space, i);
                      f(wk, i) = e(wk, i);
                    });
    } else {
      std::for_each(policy, iranges.begin(), iranges.end(),
                    [&space, &f, e](const space_size_type i) {
                      const auto& wk = getElementWorkspace(space, i);
                      f(wk, i) = e(wk, i);
                    });
    }
  }  // end of assign_scalar_impl

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      std::for_each(iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) {
                      if constexpr (use_direct_assignement) {
                        f(i) = e(i);
                      } else {
                        assign_value(f(i), e(i));
                      }
                    });
    } else {
      std::for_each(policy, iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) {
                      if constexpr (use_direct_assignement) {
                        f(i) = e(i);
                      } else {
                        assign_value(f(i), e(i));
                      }
                    });
    }
  }  // end of assign_impl

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      std::for_each(iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) {
                      const auto& wk = getElementWorkspace(space, i);
                      if constexpr (use_direct_assignement) {
                        f(wk, i) = e(wk, i);
                      } else {
                        assign_value(f(wk, i), e(wk, i));
                      }
                    });
    } else {
      std::for_each(policy, iranges.begin(), iranges.end(),
                    [&f, e](const space_size_type i) {
                      const auto& wk = getElementWorkspace(space, i);
                      if constexpr (use_direct_assignement) {
                        f(wk, i) = e(wk, i);
                      } else {
                        assign_value(f(wk, i), e(wk, i));
                      }
                    });
    }
  }  // end of assign_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void assign_sequential(FunctionType& f,
                                   const EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>) {
    using function_result_type = function_result<FunctionType>;
    using evaluator_result_type = evaluator_result<EvaluatorType>;
    constexpr auto policy = std::execution::seq;
    EvaluatorType ev = e;
    allocateWorkspace(ev);
    if constexpr ((internals::same_decay_type<function_result_type, real>)&&  //
                  (internals::same_decay_type<evaluator_result_type, real>)) {
      internals::assign_scalar_impl(policy, f, ev);
    } else {
      internals::assign_impl(policy, f, ev);
    }
  }  // end of assign_sequential

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  void assign_parallel(FunctionType& f,
                       const ExecutionPolicy policy,
                       const EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>) {
    using function_result_type = function_result<FunctionType>;
    using evaluator_result_type = evaluator_result<EvaluatorType>;
    thread_local EvaluatorType ev = e;
    allocateWorkspace(ev);
    if constexpr ((internals::same_decay_type<function_result_type, real>)&&  //
                  (internals::same_decay_type<evaluator_result_type, real>)) {
      assign_scalar_impl(policy, f, ev);
    } else {
      assign_impl(policy, f, ev);
    }
  }  // end of assign_parallel

  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  constexpr real scalar_reduce_sequential(const EvaluatorType e,
                                          const OperatorType op,
                                          const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    using Space = std::decay_t<decltype(getSpace(e))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    EvaluatorType ev = e;
    allocateWorkspace(ev);
    const auto& space = getSpace(ev);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    auto get_value = [&space, &ev](const space_size_type i) {
      if constexpr (hasElementWorkspace<Space>) {
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          const auto& wk = getElementWorkspace(space, i);
          return ev(wk, i);
        } else {
          const auto& wk = getElementWorkspace(space, i);
          return ev(wk, i)[0];
        }
      } else {
        static_cast<void>(space);
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          return ev(i);
        } else {
          return ev(i)[0];
        }
      }
    };
    return std::transform_reduce(iranges.begin(), iranges.end(), initial_value,
                                 op, get_value);
  }  // end of scalar_reduce_sequential

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  real scalar_reduce_parallel_2(const ExecutionPolicy policy,
                                const EvaluatorType ev,
                                const OperatorType op,
                                const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    using Space = std::decay_t<decltype(getSpace(ev))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(ev);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    auto get_value = [&space, ev](const space_size_type i) {
      if constexpr (hasElementWorkspace<Space>) {
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          const auto& wk = getElementWorkspace(space, i);
          return ev(wk, i);
        } else {
          const auto& wk = getElementWorkspace(space, i);
          return ev(wk, i)[0];
        }
      } else {
        static_cast<void>(space);
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          return ev(i);
        } else {
          return ev(i)[0];
        }
      }
    };
    return std::transform_reduce(policy, iranges.begin(), iranges.end(),
                                 initial_value, op, get_value);
  }

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  real scalar_reduce_parallel(const ExecutionPolicy policy,
                              const EvaluatorType& e,
                              const OperatorType op,
                              const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    thread_local EvaluatorType ev = e;
    allocateWorkspace(ev);
    return scalar_reduce_parallel_2(policy, ev, op, initial_value);
  }  // end of scalar_reduce_parallel

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType& f,
                        const ExecutionPolicy policy,
                        const EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>) {
    if (!areEquivalent(getSpace(f), getSpace(e))) {
      return ctx.registerErrorMessage("unmatched spaces");
    }
    if (getNumberOfComponents(f) != getNumberOfComponents(e)) {
      return ctx.registerErrorMessage("unmatched number of components");
    }
    //
    if (!check(ctx, e)) {
      return false;
    }
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      internals::assign_sequential(f, e);
    } else {
      internals::assign_parallel(f, policy, e);
    }
    return true;
  }  // end of assign

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  constexpr std::optional<real> scalar_reduce(AbstractErrorHandler& ctx,
                                              const ExecutionPolicy policy,
                                              const EvaluatorType e,
                                              const OperatorType op,
                                              const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    if (getNumberOfComponents(e) != 1) {
      return ctx.registerErrorMessage("non scalar evaluator");
    }
    //
    if (!check(ctx, e)) {
      return false;
    }
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      return internals::scalar_reduce_sequential(e, op, initial_value);
    } else {
      return internals::scalar_reduce_parallel(policy, e, op, initial_value);
    }
  }  // end of scalar_reduce

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType& f,
                        const EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>) {
    return assign(ctx, f, std::execution::seq, e);
  }  // end of assign

  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  constexpr std::optional<real> scalar_reduce(AbstractErrorHandler& ctx,
                                              const EvaluatorType e,
                                              const OperatorType op,
                                              const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    return scalar_reduce(ctx, std::execution::seq, e, op, initial_value);
  }  // end of scalar_reduce

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
