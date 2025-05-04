/*!
 * \file   BasicLinearQuadratureSpace.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/BasicLinearQuadratureSpace.hxx"

namespace mgis::function {

  static_assert(LinearSpaceConcept<BasicLinearQuadratureSpace<4>>);
  static_assert(QuadratureSpaceConcept<BasicLinearQuadratureSpace<4>>);
  static_assert(FunctionalSpaceConcept<BasicLinearQuadratureSpace<4>>);

}  // end of namespace mgis::function
