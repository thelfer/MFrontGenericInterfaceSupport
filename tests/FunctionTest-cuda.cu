/*!
 * \file   tests/ConstructTest.cxx
 * \brief  This file contains unit tests of `MGIS/Function`'s usage with the
 * `CUDA` programming model `make_shared` functions
 * \author Thomas Helfer
 * \date 30/06/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <thrust/universal_vector.h>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Utilities/Construct.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/CUDA/Algorithms.hxx"

template <typename Lambda>
__global__ void kernel(Lambda f) {
  f(threadIdx.x);
}

// inline constexpr auto my_negate =
//     [] __host__ __device__<typename T>(const T &v) { return -v; };

// template <typename View, typename Callable> struct MyUnaryOperationModifier {
//   MyUnaryOperationModifier(View v_, Callable c_) : v(v_), c(c_) {}
//
//   __host__ __device__ auto operator()(mgis::size_type i) const { return
//   c(v(i)); }
//
//   View v;
//   Callable c;
// };
//
// template <typename Callable> struct MyUnaryOperationModifier2 {
//   MyUnaryOperationModifier2(Callable c_) : c(c_) {}
//
//   template <typename View> __host__ __device__ auto operator()(View v) const
//   {
//     return MyUnaryOperationModifier<std::decay_t<View>, Callable>(v,
//     this->c);
//   }
//
//   Callable c;
// };
//
// template <typename Callable> auto my_unary_operation_modifier(Callable c) {
//   return MyUnaryOperationModifier2<std::decay_t<Callable>>(c);
// }
//
// inline constexpr auto my_negate =
//     [] __host__ __device__<typename View>(View v) {
//       struct NegateOperator {
//         template <typename OperandType>
//         __host__ __device__ constexpr auto
//         operator()(const OperandType &a) const noexcept
//             -> tfel::math::UnaryOperationResult<OperandType,
//                                                 tfel::math::OpNeg> //
//         requires(mgis::function::compile_time_size<
//                      tfel::math::UnaryOperationResult<OperandType,
//                                                       tfel::math::OpNeg>> !=
//                  mgis::dynamic_extent) { //
//           return -a;
//         }
//       };
//       constexpr auto nop = NegateOperator{};
//       const auto m = ::my_unary_operation_modifier(nop);
//       return m(v);
//     };

struct CUDAFunctionTest final : public tfel::tests::TestCase {
  CUDAFunctionTest()
      : tfel::tests::TestCase("MGIS/Function/CUDA", "FunctionTest") {
  }  // end of CUDAFunctionTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }  // end of execute

  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    constexpr auto layout = FunctionDataLayoutDescription{1, 1};
    auto ctx = Context{};
    auto a = thrust::universal_vector<real>{4};
    auto b = thrust::universal_vector<real>{4};
    for (mgis::size_type i = 0; i != 4; ++i) {
      b[i] = 2 * i;
    }
    auto space = BasicLinearSpace{4};
    auto a_view = FunctionView<BasicLinearSpace, layout, true>(
        space, std::span<real>(thrust::raw_pointer_cast(a.data()), 4));
    auto b_view = FunctionView<BasicLinearSpace, layout, false>(
        space, std::span<real>(thrust::raw_pointer_cast(b.data()), 4));
    auto c = CUDAExecutionConfiguration{.number_of_threads_per_block = 4};
    {
      const auto ok = assign(ctx, c, a_view, b_view | multiply_by_scalar(3));
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] - 6 * i) < eps);
      }
    }
    {
      const auto ok = assign(ctx, c, a_view, b_view | divide_by_scalar(3));
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] - static_cast<real>(2 * i) / 3) < eps);
      }
    }
  }

 private:
};

TFEL_TESTS_GENERATE_PROXY(CUDAFunctionTest, "CUDAFunctionTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FunctionTest-cuda.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
