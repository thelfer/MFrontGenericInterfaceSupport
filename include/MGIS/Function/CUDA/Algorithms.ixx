/*!
 * \file   MGIS/Function/CUDA/Algorithms.ixx
 * \brief
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

  template <typename CallableType, typename SpaceSizeType>
  __global__ void apply_kernel(CallableType f, const SpaceSizeType s) {
    const auto idx =
        static_cast<SpaceSizeType>(blockIdx.x * blockDim.x + threadIdx.x);
    if (idx >= s) {
      return;
    }
    f(idx);
  }  // end of kernel

  template <typename CallableType, typename SpaceSizeType>
  [[nodiscard]] bool call_kernel(Context& ctx,
                                 const CUDAExecutionConfiguration& c,
                                 CallableType f,
                                 const SpaceSizeType s) {
    const int nthreads = c.number_of_threads_per_block;
    const int nblocks = static_cast<int>(s + 1) / nthreads;
    if (nthreads < 1) {
      return ctx.registerErrorMessage("invalid number of threads per block");
    }
    apply_kernel<<<nblocks, nthreads>>>(f, s);
    const auto cuda_launch_error = cudaGetLastError();
    if (cuda_launch_error != cudaSuccess) {
      return ctx.registerErrorMessage(
          std::string{cudaGetErrorString(cuda_launch_error)});
    }
    const auto cuda_kernel_error = cudaDeviceSynchronize();
    if (cuda_kernel_error != cudaSuccess) {
      return ctx.registerErrorMessage(
          std::string{cudaGetErrorString(cuda_kernel_error)});
    }
    return true;
  }  // end of call_kernel

}  // end of namespace mgis::function::cuda

namespace mgis::function {

  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  bool assign(Context& ctx,
              const CUDAExecutionConfiguration& c,
              FunctionType& f,
              const EvaluatorType e)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
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
