/*!
 * \file   BasicLinearQuadratureSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"

namespace mgis::function {

  static_assert(ElementSpaceConcept<BasicLinearQuadratureSpace<4>>);
  static_assert(LinearElementSpaceConcept<BasicLinearQuadratureSpace<4>>);
  static_assert(QuadratureSpaceConcept<BasicLinearQuadratureSpace<4>>);
  static_assert(LinearQuadratureSpaceConcept<BasicLinearQuadratureSpace<4>>);

}  // end of namespace mgis::function
