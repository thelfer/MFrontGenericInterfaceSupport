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

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_scalar_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    if constexpr (internals::same_decay_type<result_type, real>) {
      for (size_type i = 0; i != ne; ++i) {
        if constexpr (internals::same_decay_type<value_type, real>) {
          f(i) = e(i);
        } else {
          f(i)[0] = e(i);
        }
      }
    } else {
      for (size_type i = 0; i != ne; ++i) {
        if constexpr (internals::same_decay_type<value_type, real>) {
          f(i) = e(i)[0];
        } else {
          f(i)[0] = e(i)[0];
        }
      }
    }
  }  // end of assign_scalar_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_scalar_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    if constexpr (internals::same_decay_type<result_type, real>) {
      for (size_type i = 0; i != ne; ++i) {
        const auto& wk = space.getElementWorkspace(i);
        if constexpr (internals::same_decay_type<value_type, real>) {
          f(wk, i) = e(wk, i);
        } else {
          f(wk, i)[0] = e(wk, i);
        }
      }
    } else {
      for (size_type i = 0; i != ne; ++i) {
        const auto& wk = space.getElementWorkspace(i);
        if constexpr (internals::same_decay_type<value_type, real>) {
          f(wk, i) = e(wk, i)[0];
        } else {
          f(wk, i)[0] = e(wk, i)[0];
        }
      }
    }
  }  // end of assign_scalar_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          !hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (use_direct_assignement) {
        f(i) = e(i);
      } else {
        auto lhs_values = f(i);
        const auto& rhs_values = e(i);
        std::copy(rhs_values.begin(), rhs_values.end(), lhs_values.begin());
      }
    }
  }  // end of assign_impl

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  void assign_impl(FunctionType& f, EvaluatorType e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto ne = getSpaceSize(space);
    for (size_type i = 0; i != ne; ++i) {
      const auto& wk = space.getElementWorkspace(i);
      if constexpr (use_direct_assignement) {
        f(wk, i) = e(wk, i);
      } else {
        auto lhs_values = f(wk, i);
        const auto& rhs_values = e(wk, i);
        std::copy(rhs_values.begin(), rhs_values.end(), lhs_values.begin());
      }
    }
  }  // end of assign_impl

}  // end of namespace mgis::function::internals

namespace mgis::function {

  //   template <size_type N, FunctionEvalutorConcept EvaluatorType>
  //   bool assign(Context& ctx, Function& f, EvaluatorType e) requires(N > 0) {
  //     checkMatchingSpaces(f, e);
  //     if (f.getNumberOfComponents() != N) {
  //       return ctx.registerErrorMessage(
  //           "assign: invalid number of components for the left hand size");
  //     }
  //     if (e.getNumberOfComponents() != N) {
  //       return ctx.registerErrorMessage(
  //           "assign: invalid number of components for the right hand size");
  //     }
  //     //
  //     if (!e.check(ctx)) {
  //       return false;
  //     }
  //     e.allocateWorkspace();
  //     //
  //     const auto qspace = getSpace(f);
  //     const auto ne = getSpaceSize(qspace);
  //     for (size_type i = 0; i != ne; ++i) {
  //       if constexpr (N == 1) {
  //         auto& v = f.getValue(i);
  //         v = e(i);
  //       } else {
  //         auto lhs_values = f.template getValues<N>(i);
  //         const auto& rhs_values = e(i);
  //         algorithm::copy<N>(rhs_values.begin(), rhs_values.end(),
  //                            lhs_values.begin());
  //       }
  //     }
  //     return true;
  //   }  // end of assign

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
    using result_type = std::invoke_result_t<EvaluatorType, size_type>;
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
    if (f.getNumberOfComponents() == 1) {
      internals::assign_scalar_impl(f, e);
    } else {
      // this if constexpr to avoid the compilation of the body
      if constexpr (!internals::same_decay_type<result_type, real>) {
        internals::assign_impl(f, e);
      }
    }
    return true;
  }  // end of assign

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
