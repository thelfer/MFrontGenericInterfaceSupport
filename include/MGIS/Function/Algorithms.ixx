/*!
 * \file   MGIS/Function/Algorithm.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_IXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_IXX

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
  void copy(const InputIterator p,
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

  inline void assign_value(auto& lhs, const auto& rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  inline void assign_value(std::span<real> lhs, const auto& rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  inline void assign_value(std::span<real> lhs,
                           const std::span<const real> rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  inline void assign_value(std::span<real> lhs, const std::span<real> rhs) {
    std::copy(rhs.begin(), rhs.end(), lhs.begin());
  }

  inline void assign_value(real& lhs, const std::span<real>& rhs) {
    lhs = rhs[0];
  }

  inline void assign_value(real& lhs, const std::span<const real>& rhs) {
    lhs = rhs[0];
  }

  inline void assign_value(real& lhs, const std::span<real, 1u>& rhs) {
    lhs = rhs[0];
  }

  inline void assign_value(real& lhs, const std::span<const real, 1u>& rhs) {
    lhs = rhs[0];
  }

  inline void assign_value(real& lhs, const std::array<real, 1u>& rhs) {
    lhs = rhs[0];
  }

  inline void assign_value(std::array<real, 1u>& lhs, const real rhs) {
    lhs[0] = rhs;
  }

  inline void assign_value(std::array<real, 1u>& lhs,
                           const std::array<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::array<real, 1u>& lhs, const std::span<real> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::array<real, 1u>& lhs, const std::span<const real> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::array<real, 1u>& lhs, const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::array<real, 1u>& lhs, const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs, const real rhs) {
    lhs[0] = rhs;
  }

  inline void assign_value(std::span<real> lhs, const real rhs) {
    lhs[0] = rhs;
  }

  inline void assign_value(std::span<real> lhs, const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real> lhs,
                           const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs, const std::span<real> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs,
                           const std::span<const real> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs,
                           const std::span<real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs,
                           const std::span<const real, 1u> rhs) {
    lhs[0] = rhs[0];
  }

  inline void assign_value(std::span<real, 1u> lhs, const std::array<real, 1u>& rhs) {
    lhs[0] = rhs[0];
  }

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_scalar_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
      f(i) = e(i);
    }
  }  // end of assign_scalar_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_scalar_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
      const auto& wk = getElementWorkspace(space, i);
      f(wk, i) = e(wk, i);
    }
  }  // end of assign_scalar_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
      if constexpr (use_direct_assignement) {
        f(i) = e(i);
      } else {
        assign_value(f(i), e(i));
      }
    }
  }  // end of assign_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
      const auto& wk = getElementWorkspace(space, i);
      if constexpr (use_direct_assignement) {
        f(wk, i) = e(wk, i);
      } else {
        assign_value(f(wk, i), e(wk, i));
      }
    }
  }  // end of assign_impl

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  bool assign(Context& ctx,
              FunctionType& f,
              EvaluatorType e)  //
      requires(((LinearElementSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>) ||
                (LinearQuadratureSpaceConcept<std::decay_t<
                     decltype(getSpace(std::declval<EvaluatorType>()))>>)) &&
               internals::same_decay_type<
                   decltype(getSpace(std::declval<FunctionType>())),
                   decltype(getSpace(std::declval<EvaluatorType>()))>) {
    using function_result_type = function_result<FunctionType>;
    using evaluator_result_type = evaluator_result<EvaluatorType>;
    if (!areEquivalent(getSpace(f), getSpace(e))) {
      return ctx.registerErrorMessage("unmatched spaces");
    }
    if (f.getNumberOfComponents() != e.getNumberOfComponents()) {
      return ctx.registerErrorMessage("unmatched number of components");
    }
    //
    if (!e.check(ctx)) {
      return false;
    }
    e.allocateWorkspace();
    if constexpr ((internals::same_decay_type<function_result_type, real>)&&  //
                  (internals::same_decay_type<evaluator_result_type, real>)) {
      internals::assign_scalar_impl(f, e);
    } else {
      internals::assign_impl(f, e);
    }
    return true;
  }  // end of assign

  template <EvaluatorConcept EvaluatorType, typename OperatorType>
  [[nodiscard]] std::optional<real> scalar_reduce(Context& ctx,
                                                  EvaluatorType e,
                                                  const OperatorType op,
                                                  const real initial_value)  //
      requires(LinearElementSpaceConcept<std::decay_t<
                   decltype(getSpace(std::declval<EvaluatorType>()))>>) {
    using Space = std::decay_t<decltype(getSpace(e))>;
    if (e.getNumberOfComponents() != 1) {
      return ctx.registerErrorMessage("non scalar evaluator");
    }
    //
    if (!e.check(ctx)) {
      return false;
    }
    e.allocateWorkspace();
    const auto& space = getSpace(e);
    const auto ne = getSpaceSize(space);
    auto r = initial_value;
    if constexpr (hasElementWorkspace<Space>) {
      for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          const auto& wk = getElementWorkspace(space, i);
          r = op(r, e(wk, i));
        } else {
          const auto& wk = getElementWorkspace(space, i);
          r = op(r, e(wk, i)[0]);
        }
      }
    } else {
      for (typename SpaceTraits<Space>::size_type i = 0; i != ne; ++i) {
        if constexpr (internals::same_decay_type<
                          evaluator_result<EvaluatorType>, real>) {
          r = op(r, e(i));
        } else {
          r = op(r, e(i)[0]);
        }
      }
    }
    return r;
  }  // end of scalar_reduce

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
