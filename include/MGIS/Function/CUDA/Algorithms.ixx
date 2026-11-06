/*!
 * \file   MGIS/Function/CUDA/Algorithms.ixx
 * \brief  Implementation of the CUDA `assign` algorithm declared in
 * `MGIS/Function/CUDA/Algorithms.hxx`, including the generic CUDA kernel and
 * the helper which launches it.
 * \author Thomas Helfer
 * \date   28/06/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_IXX
#define LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_IXX

#include "MGIS/Function/Algorithms.hxx"

namespace mgis::function::cuda {

  /*!
   * \brief CUDA kernel calling a callable object on every element of a space
   *
   * Each thread computes its own linear index from the block and thread
   * indices and calls `f` with this index. Threads whose index is greater than
   * or equal to the size of the space return immediately, which allows the
   * grid to oversubscribe the space.
   *
   * \tparam CallableType: type of the callable object, which must be callable
   * from the device
   * \tparam SpaceSizeType: integral type used to index the space
   * \param[in] f: callable object applied to every element
   * \param[in] s: size of the space
   */
  template <typename CallableType, typename SpaceSizeType>
  __global__ void apply_kernel(CallableType f, const SpaceSizeType s) {
    const auto idx =
        static_cast<SpaceSizeType>(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= s) {
      return;
    }
    f(idx);
  }  // end of kernel

  /*!
   * \brief launch `apply_kernel` on the space and wait for its completion
   *
   * The number of threads per block is taken from the configuration and the
   * number of blocks is computed from the size of the space so that every
   * index is covered by at least one thread. Both the launch and the execution
   * of the kernel are checked: the corresponding CUDA error message is
   * registered through the error handler when an error occurs.
   *
   * \tparam CallableType: type of the callable object, which must be callable
   * from the device
   * \tparam SpaceSizeType: integral type used to index the space
   * \param[in] ctx: error handler used to report the errors
   * \param[in] c: CUDA execution configuration
   * \param[in] f: callable object applied to every element of the space
   * \param[in] s: size of the space
   * \return true if the kernel was launched and executed without error
   */
  template <typename CallableType, typename SpaceSizeType>
  [[nodiscard]] bool call_kernel(AbstractErrorHandler& ctx,
                                 const CUDAExecutionConfiguration& c,
                                 CallableType f,
                                 const SpaceSizeType s) {
    const int nthreads = c.number_of_threads_per_block;
    if (nthreads < 1) {
      return ctx.registerErrorMessage("invalid number of threads per block");
    }
    const int nblocks = static_cast<int>(s + nthreads + 1) / nthreads;
    apply_kernel<<<nblocks, nthreads>>>(f, s);
    const auto cuda_launch_error = cudaGetLastError();
    if (cuda_launch_error != cudaSuccess) {
      return ctx.registerErrorMessage(cudaGetErrorString(cuda_launch_error));
    }
    const auto cuda_kernel_error = cudaDeviceSynchronize();
    if (cuda_kernel_error != cudaSuccess) {
      return ctx.registerErrorMessage(cudaGetErrorString(cuda_kernel_error));
    }
    return true;
  }  // end of call_kernel

}  // end of namespace mgis::function::cuda

namespace mgis::function {

  /*!
   * \brief implementation of the CUDA `assign` algorithm
   *
   * After checking that the spaces of the function and of the evaluator are
   * equivalent, that their numbers of components match and that the evaluator
   * is valid, a device callable lambda expression is built to assign, for a
   * given index, the value computed by the evaluator to that of the function
   * and is applied to every element of the space by `cuda::call_kernel`.
   *
   * Two orthogonal choices are made at compile time:
   * - if the function is a lightweight view, it is captured directly by the
   *   lambda expression, otherwise a view of it is computed on the host and
   *   captured,
   * - if the value returned by the function can be assigned from that of the
   *   evaluator using `operator=`, a direct assignment is used, otherwise the
   *   `internals::assign_value` helper is used, which handles the cases where
   *   the two values have different, but compatible, types.
   *
   * \note the `requires`-clause restricts the algorithm to evaluators defined
   * on a linear element space and to functions defined on the very same space
   * so that a single linear index can address both the function and the
   * evaluator.
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  bool assign(AbstractErrorHandler& ctx,
              const CUDAExecutionConfiguration& c,
              FunctionType& f,
              const EvaluatorType e)  //
      requires(
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>) {
    using Space = std::decay_t<decltype(getSpace(f))>;
    using space_size_type = typename SpaceTraits<Space>::size_type;
    using value_type = std::invoke_result_t<FunctionType, size_type>;
    using result_type = std::invoke_result_t<const EvaluatorType, size_type>;
    constexpr auto use_direct_assignement =
        requires(value_type & v1, const result_type& v2) {
      v1 = v2;
    };
    //
    const auto& space = getSpace(f);
    const auto space_size = getSpaceSize(space);
    if (!areEquivalent(space, getSpace(e))) {
      return ctx.registerErrorMessage("unmatched spaces");
    }
    if (getNumberOfComponents(f) != getNumberOfComponents(e)) {
      return ctx.registerErrorMessage("unmatched number of components");
    }
    if (!check(ctx, e)) {
      return false;
    }
    //
    if constexpr (LightweightViewConcept<FunctionType>) {
      if constexpr (use_direct_assignement) {
        auto fct =
            [f, e] MGIS_HOST_DEVICE(const space_size_type i) mutable -> void {
          f(i) = e(i);
        };
        return cuda::call_kernel(ctx, c, fct, space_size);
      } else {
        auto fct =
            [f, e] MGIS_HOST_DEVICE(const space_size_type i) mutable -> void {
          ::mgis::function::internals::assign_value(f(i), e(i));
        };
        return cuda::call_kernel(ctx, c, fct, space_size);
      }
    } else {
      auto v = view(f);
      if constexpr (use_direct_assignement) {
        auto fct =
            [v, e] MGIS_HOST_DEVICE(const space_size_type i) mutable -> void {
          v(i) = e(i);
        };
        return cuda::call_kernel(ctx, c, fct, space_size);
      } else {
        auto fct =
            [v, e] MGIS_HOST_DEVICE(const space_size_type i) mutable -> void {
          ::mgis::function::internals::assign_value(v(i), e(i));
        };
        return cuda::call_kernel(ctx, c, fct, space_size);
      }
    }
    return true;
  }  // end of assign

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_HXX */
