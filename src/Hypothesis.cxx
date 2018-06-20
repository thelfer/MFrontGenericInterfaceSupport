/*!
 * \file   Hypothesis.cxx
 * \brief
 * \author th202608
 * \date   19/06/2018
 */

#include "MFront/Behaviour/Hypothesis.hxx"
#include "MFront/Raise.hxx"

namespace mfront {

namespace behaviour {

const char *toString(const Hypothesis h) {
  if (h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) {
    return "AxisymmetricalGeneralisedPlaneStrain";
  } else if (h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS) {
    return "AxisymmetricalGeneralisedPlaneStress";
  } else if (h == Hypothesis::AXISYMMETRICAL) {
    return "Axisymmetrical";
  } else if (h == Hypothesis::PLANESTRESS) {
    return "PlaneStress";
  } else if (h == Hypothesis::PLANESTRAIN) {
    return "PlaneStrain";
  } else if (h == Hypothesis::GENERALISEDPLANESTRAIN) {
    return "GeneralisedPlaneStrain";
  } else if (h == Hypothesis::TRIDIMENSIONAL) {
    return "Tridimensional";
  }
  raise("toString : unsupported modelling hypothesis");
} // end of toString

} // end of namespace behaviour

} // end of namespace mfront
