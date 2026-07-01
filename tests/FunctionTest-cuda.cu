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
#include "MGIS/Function/TFEL/Tensors.hxx"
#include "MGIS/Function/TFEL/Mechanics.hxx"
#include "MGIS/Function/CUDA/Algorithms.hxx"

template <typename Lambda>
__global__ void kernel(Lambda f) {
  f(threadIdx.x);
}

struct CUDAFunctionTest final : public tfel::tests::TestCase {
  CUDAFunctionTest()
      : tfel::tests::TestCase("MGIS/Function/CUDA", "FunctionTest") {
  }  // end of CUDAFunctionTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }  // end of execute

  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    constexpr auto layout = FunctionDataLayoutDescription{1, 1};
    constexpr auto layout2 =
        FunctionDataLayoutDescription{.data_size = 2, .data_stride = 2};
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
    auto conf = CUDAExecutionConfiguration{.number_of_threads_per_block = 4};
    {
      const auto ok = assign(ctx, conf, a_view, b_view | multiply_by_scalar(3));
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] - 6 * i) < eps);
      }
    }
    {
      const auto ok = assign(ctx, conf, a_view, b_view | divide_by_scalar(3));
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] - static_cast<real>(2 * i) / 3) < eps);
      }
    }
    {
      const auto ok = assign(ctx, conf, a_view, b_view | negate);
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] + 2 * i) < eps);
      }
    }
    //
    for (mgis::size_type i = 0; i != 4; ++i) {
      b[i] = (i % 2 == 0)  //
                 ? 2 * static_cast<real>(i)
                 : -3 * static_cast<real>(i);
    }
    {
      const auto ok = assign(ctx, conf, a_view, b_view | absolute_value);
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        if (i % 2 == 0) {
          TFEL_TESTS_ASSERT(std::abs(a[i] - 2 * i) < eps);
        } else {
          TFEL_TESTS_ASSERT(std::abs(a[i] - 3 * i) < eps);
        }
      }
    }
    //
    auto c = thrust::universal_vector<real>{4 * 2};
    for (mgis::size_type i = 0; i != 4; ++i) {
      c[i * 2] = 2 * i;
      c[i * 2 + 1] = -3 * static_cast<real>(i);
    }
    auto c_view = FunctionView<BasicLinearSpace, layout2, false>(
        space, std::span<const real>(thrust::raw_pointer_cast(c.data()), 8));
    {
      const auto ok = assign(ctx, conf, a_view, c_view | maximum_component);
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] - 2 * i) < eps);
      }
    }
    {
      const auto ok = assign(ctx, conf, a_view, c_view | minimum_component);
      TFEL_TESTS_ASSERT(ok);
      for (mgis::size_type i = 0; i != 4; ++i) {
        TFEL_TESTS_ASSERT(std::abs(a[i] + 3 * i) < eps);
      }
    }
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto conf = CUDAExecutionConfiguration{.number_of_threads_per_block = 4};
    auto ctx = Context{};
    auto space = BasicLinearSpace{2};
    auto values = thrust::universal_vector<real>{1, 2, 3, 4};
    const auto f =
        FunctionEvaluator<BasicLinearSpace>(
            space, std::span<real>(thrust::raw_pointer_cast(values.data()), 4),
            2) |
        as_tvector<2>;
    TFEL_TESTS_ASSERT(f.check(ctx));
    auto values2 = thrust::universal_vector<real>{1, -2, 3, -4};
    const auto f2 =
        FunctionEvaluator<BasicLinearSpace>(
            space, std::span<real>(thrust::raw_pointer_cast(values2.data()), 4),
            2) |
        as_tvector<2>;
    TFEL_TESTS_ASSERT(f2.check(ctx));
    auto results = thrust::universal_vector<real>{4};
    auto rview =
        FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, true>(
            space, std::span<real>(thrust::raw_pointer_cast(results.data()), 4),
            2);
    auto r = rview | as_tvector<2>;
    auto results2 = thrust::universal_vector<real>{2};
    auto r2 = FunctionView<BasicLinearSpace,
                           FunctionDataLayoutDescription{1, 1}, true>(
        space, std::span<real>(thrust::raw_pointer_cast(results.data()), 4));
    TFEL_TESTS_ASSERT(r.check(ctx));
    auto reset = [&results, &results2] {
      std::fill(results.begin(), results.end(), 0.);
      std::fill(results2.begin(), results2.end(), 0.);
    };
    //
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | negate);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] + 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] + 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] + 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] + 4) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | multiply_by_scalar(4));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 4) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] - 8) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 12) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] - 16) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | divide_by_scalar(3));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(3 * r(0)[0] - 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(3 * r(0)[1] - 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(3 * r(1)[0] - 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(3 * r(1)[1] - 4) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, add(f, f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] - 0) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 6) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] - 0) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | add(f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] - 0) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 6) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] - 0) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | substract(f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 0) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] - 4) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 0) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] - 8) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | mean_value(f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] - 0) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] - 0) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r2, f | inner_product(f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r2(0) + 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r2(1) + 7) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r2, inner_product(f, f2));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r2(0) + 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r2(1) + 7) < eps);
    }
    {
      reset();
      const auto ok = assign(ctx, conf, r, f | divide(f | inner_product(f2)));
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(3 * r(0)[0] + 1) < eps);
      TFEL_TESTS_ASSERT(std::abs(3 * r(0)[1] + 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(7 * r(1)[0] + 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(7 * r(1)[1] + 4) < eps);
    }
    {
      reset();
      const auto op = binary_operation(
          [] MGIS_HOST_DEVICE(const tfel::math::tvector<2u, real>& v1,
                              const tfel::math::tvector<2u, real>& v2) {
            return eval(v1 + 2 * v2);
          },
          f, f2);
      const auto ok = assign(ctx, conf, r, op);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] + 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 9) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] + 4) < eps);
    }
    {
      reset();
      const auto op =
          f | binary_operation(
                  [] MGIS_HOST_DEVICE(const tfel::math::tvector<2u, real>& v1,
                                      const tfel::math::tvector<2u, real>& v2) {
                    return eval(v1 + 2 * v2);
                  },
                  f2);
      const auto ok = assign(ctx, conf, r, op);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] + 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 9) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] + 4) < eps);
    }
    {
      reset();
      const auto op = binary_operation(
          [] MGIS_HOST_DEVICE(const tfel::math::tvector<2u, real>& v1,
                              const tfel::math::tvector<2u, real>& v2) {
            return eval(v1 + 2 * v2);
          },
          f, f2);
      const auto ok = assign(ctx, conf, r, op);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(r(0)[0] - 3) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(0)[1] + 2) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[0] - 9) < eps);
      TFEL_TESTS_ASSERT(std::abs(r(1)[1] + 4) < eps);
    }
  }
};

struct CUDAMechanicalEvaluatorsTest final : public tfel::tests::TestCase {
  CUDAMechanicalEvaluatorsTest()
      : tfel::tests::TestCase("MGIS/Function/CUDA", "MechanicalEvaluatorsTests") {
  }  // end of CUDAMechanicalEvaluatorsTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    return this->result;
  }

 private:
  void test1() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    const auto conf =
        CUDAExecutionConfiguration{.number_of_threads_per_block = 32};
    auto ctx = Context{};
    auto space = BasicLinearSpace(1);
    auto values = thrust::universal_vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(
        space, std::span<real>(thrust::raw_pointer_cast(values.data()), 3), 3);
    TFEL_TESTS_CHECK_EQUAL(f.getNumberOfComponents(), 3);
    TFEL_TESTS_CHECK_EQUAL(f.getDataStride(), 3);
    const auto s = view<3>(f) | as_stensor<1>;
    TFEL_TESTS_ASSERT(s.check(ctx));
    //
    auto results = thrust::universal_vector<real>{1};
    auto r =
        FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, true>(
            space, std::span<real>(thrust::raw_pointer_cast(results.data()), 1),
            1);
    {
      const auto e = vmis(s);
      const auto ok = assign(ctx, conf, r, e);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(results[0] - std::sqrt(3)) < eps);
    }
    {
      const auto e = hydrostatic_stress(s);
      const auto ok = assign(ctx, conf, r, e);
      TFEL_TESTS_ASSERT(ok);
      TFEL_TESTS_ASSERT(std::abs(results[0] - 2) < eps);
    }
  }
  void test2() {
    using namespace mgis;
    using namespace mgis::function;
    constexpr auto eps = real{1e-14};
    auto ctx = Context{};
    const auto conf =
        CUDAExecutionConfiguration{.number_of_threads_per_block = 32};
    auto space = BasicLinearSpace(1);
    const auto values = thrust::universal_vector<real>{1, 2, 3};
    const auto f = FunctionEvaluator<BasicLinearSpace>(
        space,
        std::span<const real>(thrust::raw_pointer_cast(values.data()), 3), 3);
    auto seq_results = thrust::universal_vector<real>{1};
    auto seq =
        FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, true>(
            space,
            std::span<real>(thrust::raw_pointer_cast(seq_results.data()), 1),
            1);
    TFEL_TESTS_ASSERT(seq.isScalar());
    const auto op = f | as_stensor<1> | vmis;
    const auto ok = assign(ctx, conf, seq, op);
    TFEL_TESTS_ASSERT(ok);
    TFEL_TESTS_ASSERT(std::abs(seq_results[0] - std::sqrt(3)) < eps);
  }
};

TFEL_TESTS_GENERATE_PROXY(CUDAFunctionTest, "CUDAFunctionTest");
TFEL_TESTS_GENERATE_PROXY(CUDAMechanicalEvaluatorsTest, "CUDAMechanicalEvaluatorsTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("FunctionTest-cuda.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
