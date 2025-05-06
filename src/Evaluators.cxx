/*!
 * \file   Evaluators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Evaluators.hxx"

namespace mgis::function {

  static_assert(std::same_as<
                decltype(std::declval<FixedSizeEvaluator<BasicLinearSpace, 9>>()
                             .getSpace()
                             .size()),
                size_type>);
  static_assert(std::same_as<BasicLinearSpace::size_type, size_type>);

}  // end of namespace mgis::function