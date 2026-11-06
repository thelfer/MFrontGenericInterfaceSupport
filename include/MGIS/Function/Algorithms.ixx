/*!
 * \file   MGIS/Function/Algorithms.ixx
 * \brief Implementation of the algorithms which can be applied to the functions
 * and the evaluators of the `MGIS/Function` library.
 * \author Thomas Helfer
 * \date   29/04/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_IXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_IXX

#include <ranges>
#include <numeric>
#include <iterator>
#include <algorithm>

namespace mgis::function::algorithm {

  /*!
   * \brief copies `N` values from a range of input iterators to a range of
   * output iterators
   *
   * When both iterators are random access iterators and the number of values is
   * small, the copy is unrolled at compile time, which generates a more
   * efficient code than the call to `std::copy` used otherwise.
   *
   * \tparam N: number of values to be copied, which must be greater than zero
   * \tparam InputIterator: input iterator type
   * \tparam OutputIterator: output iterator type
   * \param[in] p: iterator to the first value to be copied
   * \param[in] pe: iterator past the last value to be copied
   * \param[in] po: iterator to the first value assigned
   * \note `std::copy` is used if more than nine values are copied or if one of
   * the iterators is not a random access iterator
   */
  template <size_type N, typename InputIterator, typename OutputIterator>
  constexpr void copy(const InputIterator p,
                      const InputIterator pe,
                      const OutputIterator po) requires(N > 0) {
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
        *po = *p;
        copy<N - 1>(std::next(p), pe, std::next(po));
      }
    } else {
      std::copy(p, pe, po);
    }
  }  // end of copy

}  // end of namespace mgis::function::algorithm

namespace mgis::function::internals {

  /*!
   * \brief implementation of the most generic `assign_value` function, which
   * iterates over the values of the right hand side
   *
   * The other overloads, declared in `MGIS/Function/Algorithms.hxx`, copy the
   * values one by one, those connecting a scalar value and a container of a
   * single value using only the first value of this container.
   */
  constexpr void assign_value(std::span<real> lhs, const auto& rhs) {
    std::size_t i = 0;
    for (const auto& v : rhs) {
      lhs[i++] = v;
    }
  }

  constexpr void assign_value(std::span<real> lhs,
                              const std::span<const real> rhs) {
    for (std::size_t i = 0; i != rhs.size(); ++i) {
      lhs[i] = rhs[i];
    }
  }

  constexpr void assign_value(std::span<real> lhs, const std::span<real> rhs) {
    for (std::size_t i = 0; i != rhs.size(); ++i) {
      lhs[i] = rhs[i];
    }
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

  /*!
   * \brief assigns to a function the scalar values computed by an evaluator,
   * the evaluation being sequential
   *
   * This overload is used when the values returned by both objects are scalar
   * values and when the space does not require any workspace to access its
   * elements: the function and the evaluator are called with the index of each
   * element of the space.
   *
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   */
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

  /*!
   * \brief assigns to a function the scalar values computed by an evaluator
   * using the given execution policy
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`, which keeps this function `constexpr`
   * in that case
   * \note a view of the function is built before the loop if the function is
   * not its own view
   */
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

  /*!
   * \brief assigns to a function the scalar values computed by an evaluator,
   * the evaluation being sequential
   *
   * This overload is used when the values returned by both objects are scalar
   * values and when the space requires a workspace to access its elements: the
   * workspace of an element is retrieved before calling the function and the
   * evaluator with this workspace and the index of the element. The function
   * assigned and the evaluator providing the values assigned are the same
   * arguments as for the other overload.
   */
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

  /*!
   * \brief assigns to a function the scalar values computed by an evaluator
   * using the given execution policy
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`, which keeps this function `constexpr`
   * in that case
   * \note the workspace of an element is retrieved within the loop, as it may
   * depend on the index of the element
   */
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

  /*!
   * \brief assigns to a function the values computed by an evaluator, the
   * evaluation being sequential
   *
   * This overload is used when the values returned by the function or the
   * evaluator are not scalar values and when the space does not require any
   * workspace to access its elements.
   *
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the values returned by the evaluator are assigned directly to those
   * returned by the function if the corresponding types are compatible, using
   * one of the `assign_value` helpers otherwise
   */
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

  /*!
   * \brief assigns to a function the values computed by an evaluator using the
   * given execution policy
   *
   * This overload is used when the values returned by the function or the
   * evaluator are not scalar values and when the space does not require any
   * workspace to access its elements.
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`, which keeps this function `constexpr`
   * in that case
   * \note a view of the function is built before the loop if the function is
   * not its own view
   */
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

  /*!
   * \brief assigns to a function the values computed by an evaluator, the
   * evaluation being sequential
   *
   * This overload is used when the values returned by the function or the
   * evaluator are not scalar values and when the space requires a workspace to
   * access its elements. The function assigned and the evaluator providing the
   * values assigned are the same arguments as for the other overload.
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void
  assign_sequential_impl(FunctionType& f, const EvaluatorType& e) requires(
      (LinearElementSpaceConcept<std::decay_t<decltype(getSpace(f))>>)&&(
          hasElementWorkspace<std::decay_t<decltype(getSpace(f))>>)) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using value_type =
        std::invoke_result_t<FunctionType&, const element_workspace<Space>&,
                             size_type>;
    using result_type =
        std::invoke_result_t<const EvaluatorType&,
                             const element_workspace<Space>&, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    const auto& space = getSpace(f);
    const auto iranges =
        std::views::iota(space_size_type{}, getSpaceSize(space));
    std::for_each(iranges.begin(), iranges.end(),
                  [&space, &f, e](const space_size_type i) {
                    const auto& wk = getElementWorkspace(space, i);
                    if constexpr (use_direct_assignement) {
                      f(wk, i) = e(wk, i);
                    } else {
                      assign_value(f(wk, i), e(wk, i));
                    }
                  });
  }  // end of assign_sequential_impl

#ifdef MGIS_HAS_STL_PARALLEL_ALGORITHMS

  /*!
   * \brief assigns to a function the values computed by an evaluator using the
   * given execution policy
   *
   * This overload is used when the values returned by the function or the
   * evaluator are not scalar values and when the space requires a workspace to
   * access its elements.
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`, which keeps this function `constexpr`
   * in that case
   * \note a view of the function is built before the loop if the function is
   * not its own view
   */
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
      using value_type =
          std::invoke_result_t<FunctionType&, const element_workspace<Space>&,
                               size_type>;
      using result_type =
          std::invoke_result_t<const EvaluatorType&,
                               const element_workspace<Space>&, size_type>;
      constexpr auto use_direct_assignement =
          requires(value_type & v1, const result_type& v2) {
        v1 = v2;
      };
      const auto& space = getSpace(f);
      const auto iranges =
          std::views::iota(space_size_type{}, getSpaceSize(space));
      if constexpr (LightweightViewConcept<FunctionType>) {
        std::for_each(policy, iranges.begin(), iranges.end(),
                      [&space, f, e](const space_size_type i) mutable {
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
                      [&space, v, e](const space_size_type i) mutable {
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

  /*!
   * \brief assigns to a function the values computed by an evaluator, the
   * evaluation being sequential
   *
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in,out] f: function assigned
   * \param[in] e: evaluator providing the values assigned
   * \note the scalar implementation is used when both the function and the
   * evaluator return scalar values, the general implementation otherwise
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  constexpr void assign_sequential(FunctionType& f,
                                   const EvaluatorType e)  //
      requires(
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) &&
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

  /*!
   * \brief assigns to a function the values computed by an evaluator using the
   * given execution policy
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam FunctionType: function type
   * \tparam EvaluatorType: evaluator type
   * \param[in,out] f: function assigned
   * \param[in] policy: execution policy used to perform the assignment
   * \param[in] e: evaluator providing the values assigned
   * \note the scalar implementation is used when both the function and the
   * evaluator return scalar values, the general implementation otherwise
   */
  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            typename FunctionType,
            EvaluatorConcept EvaluatorType>
  void assign_parallel(FunctionType& f,
                       const ExecutionPolicy policy,
                       const EvaluatorType e)  //
      requires((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>)&&  //
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

  /*!
   * \brief reduces the values computed by an evaluator to a single scalar
   * value, the evaluation being sequential
   *
   * The reduction is performed by `std::transform_reduce`: the values computed
   * by the evaluator are combined two by two with the reduction operator,
   * starting from the initial value.
   *
   * \tparam EvaluatorType: evaluator type
   * \tparam OperatorType: reduction operator type
   * \param[in] e: evaluator reduced
   * \param[in] op: reduction operator
   * \param[in] initial_value: value the reduction starts from
   * \return the reduced value
   * \note the first component of the values returned by the evaluator is used
   * when those values are not scalar values
   */
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

  /*!
   * \brief reduces the values computed by an evaluator to a single scalar
   * value using the given execution policy
   *
   * The reduction is performed by `std::transform_reduce`, which is given the
   * execution policy so that the evaluation can be parallelized.
   *
   * \tparam ExecutionPolicy: execution policy type
   * \tparam EvaluatorType: evaluator type
   * \tparam OperatorType: reduction operator type
   * \param[in] policy: execution policy used to perform the reduction
   * \param[in] e: evaluator reduced
   * \param[in] op: reduction operator
   * \param[in] initial_value: value the reduction starts from
   * \return the reduced value
   * \note the first component of the values returned by the evaluator is used
   * when those values are not scalar values
   */
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

  /*!
   * \brief implementation of the `assign` algorithm using an execution policy
   *
   * The error messages "unmatched spaces" and "unmatched number of
   * components" are registered in the error handler when the spaces of the
   * function and the evaluator are not equivalent or when the number of
   * components of the values they return differs. The evaluator is then
   * checked before the values are assigned.
   *
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`
   */
  template <ExecutionPolicyConceptConcept ExecutionPolicy,
            ViewableFunctionArgumentConcept FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType&& f,
                        const ExecutionPolicy policy,
                        const EvaluatorType e)  //
      requires(
          (internals::isEvaluatorAssignableToFunction<
              EvaluatorType,
              std::decay_t<FunctionType>>)&&  //
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>)) {
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

  /*!
   * \brief implementation of the sequential `assign` algorithm
   *
   * The error messages "unmatched spaces" and "unmatched number of
   * components" are registered in the error handler when the spaces of the
   * function and the evaluator are not equivalent or when the number of
   * components of the values they return differs. The evaluator is then
   * checked before the values are assigned.
   */
  template <ViewableFunctionArgumentConcept FunctionType,
            EvaluatorConcept EvaluatorType>
  constexpr bool assign(AbstractErrorHandler& ctx,
                        FunctionType&& f,
                        const EvaluatorType e)  //
      requires(
          (internals::isEvaluatorAssignableToFunction<
              EvaluatorType,
              std::decay_t<FunctionType>>)&&  //
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>)) {
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

  /*!
   * \brief implementation of the `scalar_reduce` algorithm using an execution
   * policy
   *
   * The error message "non scalar evaluator" is registered in the error
   * handler when the values returned by the evaluator have more than one
   * component. The evaluator is then checked before being reduced.
   *
   * \note the sequential implementation is called when the policy is
   * `std::execution::sequenced_policy`
   */
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
      return std::nullopt;
    }
    if constexpr (std::same_as<ExecutionPolicy,
                               std::execution::sequenced_policy>) {
      return internals::sequential_scalar_reduce(e, op, initial_value);
    } else {
      return internals::scalar_reduce_parallel(policy, e, op, initial_value);
    }
  }  // end of scalar_reduce

#endif /* MGIS_HAS_STL_PARALLEL_ALGORITHMS */

  /*!
   * \brief implementation of the sequential `scalar_reduce` algorithm
   *
   * The error message "non scalar evaluator" is registered in the error
   * handler when the values returned by the evaluator have more than one
   * component. The evaluator is then checked before being reduced.
   */
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
      return std::nullopt;
    }
    return internals::sequential_scalar_reduce(e, op, initial_value);
  }  // end of scalar_reduce

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
