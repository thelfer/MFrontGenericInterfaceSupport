/*!
 * \file   FixedSizeView.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/10/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/FixedSizeView.hxx"

namespace mgis::function {

  static_assert(FunctionConcept<FixedSizeView<Function<BasicLinearSpace>, 1>>);

}  // end of namespace mgis::function
