/*!
 * \file   Hypothesis.cxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
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
