/*!
 * \file   Evaluators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/FixedSizeEvaluator.hxx"

namespace mgis::function {

  static_assert(std::same_as<
                decltype(std::declval<FixedSizeEvaluator<BasicLinearSpace, 9>>()
                             .getSpace()
                             .size()),
                size_type>);
  static_assert(std::same_as<BasicLinearSpace::size_type, size_type>);

  static_assert(
      std::same_as<evaluator_result<FixedSizeEvaluator<BasicLinearSpace, 9>>,
                   std::span<const real, 9>>);
  static_assert(internals::CompileTimeSize<real>::value == 1);
  static_assert(internals::CompileTimeSize<std::span<real, 9>>::value == 9);
  static_assert(internals::CompileTimeSize<std::span<const real, 9>>::value ==
                9);
  static_assert(compile_time_size<std::span<const real>> == dynamic_extent);

  static_assert(number_of_components<FixedSizeEvaluator<BasicLinearSpace, 9>> ==
                9);

}  // end of namespace mgis::function
