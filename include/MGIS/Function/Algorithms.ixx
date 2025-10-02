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

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void
  assign_sequential_scalar_impl(FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    std::for_each(iranges.begin(), iranges.end(),
                  [&f, e](const space_size_type i) { f(i) = e(i); });
  }  // end of assign_sequential_scalar_impl

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_scalar_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    // this test is made for assign_scalar_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      assign_sequential_scalar_impl(f, e);
    } else {
      using Space = std::decay_t<decltype(getSpace(f))>;
      using space_size_type = typename SpaceTraits<Space>::size_type;
      const auto& space = getSpace(f);
      if constexpr (LightweightViewConcept<FunctionType>) {
        const auto iranges =
            std::views::iota(space_size_type{}, getSpaceSize(space));
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [f, e](const space_size_type i) mutable { f(i) = e(i); });
      } else {
        auto v = view(f);
        const auto iranges =
            std::views::iota(space_size_type{}, getSpaceSize(space));
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [v, e](const space_size_type i) mutable { v(i) = e(i); });
      }
    }
  }  // end of assign_scalar_impl

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void
  assign_sequential_scalar_impl(FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    std::for_each(iranges.begin(), iranges.end(),
                  [&space, &f, e](const space_size_type i) {
                    const auto& wk = getElementWorkspace(space, i);
                    f(wk, i) = e(wk, i);
                  });
  }  // end of assign_sequential_scalar_impl

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_scalar_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      assign_sequential_scalar_impl(f, e);
    } else {
      using Space = std::decay_t<decltype(getSpace(f))>;
      using space_size_type = typename SpaceTraits<Space>::size_type;
      const auto& space = getSpace(f);
      const auto iranges =
          std::views::iota(space_size_type{}, getSpaceSize(space));
      if constexpr (LightweightViewConcept<FunctionType>) {
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [&space, f, e](const space_size_type i) mutable {
                        const auto& wk = getElementWorkspace(space, i);
                        f(wk, i) = e(wk, i);
                      });
      } else {
        auto v = view(f);
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [&space, v, e](const space_size_type i) mutable {
                        const auto& wk = getElementWorkspace(space, i);
                        v(wk, i) = e(wk, i);
                      });
      }
    }
  }  // end of assign_scalar_impl

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void
  assign_sequential_impl(FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using result_type = std::invoke_result_t<const EvaluatorType, size_type>;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    std::for_each(iranges.begin(), iranges.end(),
                  [&f, e](const space_size_type i) {
                    if constexpr (use_direct_assignement) {
                      f(i) = e(i);
                    } else {
                      assign_value(f(i), e(i));
                    }
                  });
  }  // end of assign_sequential_impl

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      assign_sequential_impl(f, e);
    } else {
      using Space = std::decay_t<decltype(getSpace(f))>;
      using space_size_type = typename SpaceTraits<Space>::size_type;
      using value_type = std::invoke_result_t<FunctionType, size_type>;
      using result_type = std::invoke_result_t<const EvaluatorType, size_type>;
      constexpr auto use_direct_assignement =
          requires(value_type & v1, const result_type& v2) {
        v1 = v2;
      };
      const auto& space = getSpace(f);
      const auto iranges =
          std::views::iota(space_size_type{}, getSpaceSize(space));
      if constexpr (LightweightViewConcept<FunctionType>) {
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [f, e](const space_size_type i) mutable {
                        if constexpr (use_direct_assignement) {
                          f(i) = e(i);
                        } else {
                          assign_value(f(i), e(i));
                        }
                      });
      } else {
        auto v = view(f);
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [v, e](const space_size_type i) mutable {
                        if constexpr (use_direct_assignement) {
                          v(i) = e(i);
                        } else {
                          assign_value(v(i), e(i));
                        }
                      });
      }
    }
  }  // end of assign_impl

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void
  assign_sequential_impl(FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<const EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    std::for_each(iranges.begin(), iranges.end(),
                  [&f, e](const space_size_type i) {
                    const auto& wk = getElementWorkspace(space, i);
                    if constexpr (use_direct_assignement) {
                      f(wk, i) = e(wk, i);
                    } else {
                      assign_value(f(wk, i), e(wk, i));
                    }
                  });
  }  // end of assign_sequential_impl

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr void
  assign_impl(const ExecutionPolicy policy, FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    // this test is made for assign_impl be constexpr for non-parallel
    // evaluation
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      assign_sequential_impl(f, e);
    } else {
      using Space = std::decay_t<decltype(getSpace(f))>;
      using space_size_type = typename SpaceTraits<Space>::size_type;
      using value_type = std::invoke_result_t<FunctionType, size_type>;
      using result_type = std::invoke_result_t<const EvaluatorType, size_type>;
      constexpr auto use_direct_assignement =
          requires(value_type & v1, const result_type& v2) {
        v1 = v2;
      };
      const auto& space = getSpace(f);
      const auto iranges =
          std::views::iota(space_size_type{}, getSpaceSize(space));
      if constexpr (LightweightViewConcept<FunctionType>) {
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [f, e](const space_size_type i) mutable {
                        const auto& wk = getElementWorkspace(space, i);
                        if constexpr (use_direct_assignement) {
                          f(wk, i) = e(wk, i);
                        } else {
                          assign_value(f(wk, i), e(wk, i));
                        }
                      });
      } else {
        auto v = view(f);
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [v, e](const space_size_type i) mutable {
                        const auto& wk = getElementWorkspace(space, i);
                        if constexpr (use_direct_assignement) {
                          v(wk, i) = e(wk, i);
                        } else {
                          assign_value(v(wk, i), e(wk, i));
                        }
                      });
      }
    }
  }  // end of assign_impl

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void assign_sequential(FunctionType& f,
                                   const EvaluatorType e)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          same_decay_type<function_space<FunctionType>,
                          evaluator_space<EvaluatorType>>) {
    using function_result_type = function_result<FunctionType>;
    using evaluator_result_type = evaluator_result<EvaluatorType>;
    if constexpr ((same_decay_type<function_result_type, real>)&&  //
                  (same_decay_type<evaluator_result_type, real>)) {
      assign_sequential_scalar_impl(f, e);
    } else {
      assign_sequential_impl(f, e);
    }
  }  // end of assign_sequential

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  void assign_parallel(FunctionType& f,
                       const ExecutionPolicy policy,
                       const EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<

                     evaluator_space<EvaluatorType>>>)) &&
               same_decay_type<function_space<FunctionType>,
                               evaluator_space<EvaluatorType>>) {
    using function_result_type = function_result<FunctionType>;
    using evaluator_result_type = evaluator_result<EvaluatorType>;
    if constexpr ((same_decay_type<function_result_type, real>)&&  //
                  (same_decay_type<evaluator_result_type, real>)) {
      assign_scalar_impl(policy, f, e);
    } else {
      assign_impl(policy, f, e);
    }
  }  // end of assign_parallel

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  constexpr real sequential_scalar_reduce(const EvaluatorType e,
                                          const OperatorType op,
                                          const real initial_value)  //
      requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) {
    using Space = std::decay_t<decltype(getSpace(e))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(e);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    // warning, taking a reference to space would probably lead to
    // a segfault when offloading to GPUs
    auto get_value = [&space, e](const space_size_type i) {
      if constexpr (hasElementWorkspace<Space>) {
        if constexpr (same_decay_type<evaluator_result<EvaluatorType>, real>) {
          const auto& wk = getElementWorkspace(space, i);
          return e(wk, i);
        } else {
          const auto& wk = getElementWorkspace(space, i);
          return e(wk, i)[0];
        }
      } else {
        static_cast<void>(space);
        if constexpr (same_decay_type<evaluator_result<EvaluatorType>, real>) {
          return e(i);
        } else {
          return e(i)[0];
        }
      }
    };
    return std::transform_reduce(iranges.begin(), iranges.end(), initial_value,
                                 op, get_value);
  }  // end of sequential_scalar_reduce

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  real scalar_reduce_parallel(const ExecutionPolicy policy,
                              const EvaluatorType& e,
                              const OperatorType op,
                              const real initial_value)  //
      requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) {
    using Space = std::decay_t<decltype(getSpace(e))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    const auto& space = getSpace(e);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    auto get_value = [&space, e](const space_size_type i) {
      if constexpr (hasElementWorkspace<Space>) {
        if constexpr (same_decay_type<evaluator_result<EvaluatorType>, real>) {
          const auto& wk = getElementWorkspace(space, i);
          return e(wk, i);
        } else {
          const auto& wk = getElementWorkspace(space, i);
          return e(wk, i)[0];
        }
      } else {
        static_cast<void>(space);
        if constexpr (same_decay_type<evaluator_result<EvaluatorType>, real>) {
          return e(i);
        } else {
          return e(i)[0];
        }
      }
    };
    return std::transform_reduce(policy, iranges.begin(), iranges.end(),
                                 initial_value, op, get_value);
  }  // end of scalar_reduce_parallel

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

}  // end of namespace mgis::function::internals

namespace mgis::function {

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType& f,
                        const ExecutionPolicy policy,
                        const EvaluatorType e)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>) {
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

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType& f,
                        const EvaluatorType e)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>) {
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
    internals::assign_sequential(f, e);
    return true;
  }  // end of assign

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            EvaluatorConcept EvaluatorType,
            typename OperatorType>
  constexpr std::optional<real> scalar_reduce(AbstractErrorHandler& ctx,
                                              const ExecutionPolicy policy,
                                              const EvaluatorType e,
                                              const OperatorType op,
                                              const real initial_value)  //
      requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) {
    if (getNumberOfComponents(e) != 1) {
      return ctx.registerErrorMessage("non scalar evaluator");
    }
    //
    if (!check(ctx, e)) {
      return false;
    }
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      return internals::sequential_scalar_reduce(e, op, initial_value);
    } else {
      return internals::scalar_reduce_parallel(policy, e, op, initial_value);
    }
  }  // end of scalar_reduce

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  constexpr std::optional<real> scalar_reduce(AbstractErrorHandler& ctx,
                                              const EvaluatorType e,
                                              const OperatorType op,
                                              const real initial_value)  //
      requires(LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) {
    if (getNumberOfComponents(e) != 1) {
      return ctx.registerErrorMessage("non scalar evaluator");
    }
    //
    if (!check(ctx, e)) {
      return false;
    }
    return internals::sequential_scalar_reduce(e, op, initial_value);
  }  // end of scalar_reduce

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
