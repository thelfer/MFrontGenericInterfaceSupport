/*!
 * \file   BasicLinearSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"

namespace mgis::function {

  static_assert(ElementSpaceConcept<BasicLinearSpace>);
  static_assert(LinearElementSpaceConcept<BasicLinearSpace>);

}  // end of namespace mgis::function
