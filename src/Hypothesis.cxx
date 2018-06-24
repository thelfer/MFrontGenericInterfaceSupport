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

#include <cstring>
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Raise.hxx"

namespace mgis {

  namespace behaviour {

    Hypothesis fromString(const std::string& h) {
      if (h == "AxisymmetricalGeneralisedPlaneStrain") {
        return Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
      } else if (h == "AxisymmetricalGeneralisedPlaneStress") {
        return Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS;
      } else if (h == "Axisymmetrical") {
        return Hypothesis::AXISYMMETRICAL;
      } else if (h == "PlaneStress") {
        return Hypothesis::PLANESTRESS;
      } else if (h == "PlaneStrain") {
        return Hypothesis::PLANESTRAIN;
      } else if (h == "GeneralisedPlaneStrain") {
        return Hypothesis::GENERALISEDPLANESTRAIN;
      } else if (h == "Tridimensional") {
        return Hypothesis::TRIDIMENSIONAL;
      }
      raise("toString : unsupported modelling hypothesis");
    }  // end of fromString

    Hypothesis fromString(const char* const h) {
      if (::strcmp(h, "AxisymmetricalGeneralisedPlaneStrain") == 0) {
        return Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
      } else if (::strcmp(h, "AxisymmetricalGeneralisedPlaneStress") == 0) {
        return Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS;
      } else if (::strcmp(h, "Axisymmetrical") == 0) {
        return Hypothesis::AXISYMMETRICAL;
      } else if (::strcmp(h, "PlaneStress") == 0) {
        return Hypothesis::PLANESTRESS;
      } else if (::strcmp(h, "PlaneStrain") == 0) {
        return Hypothesis::PLANESTRAIN;
      } else if (::strcmp(h, "GeneralisedPlaneStrain") == 0) {
        return Hypothesis::GENERALISEDPLANESTRAIN;
      } else if (::strcmp(h, "Tridimensional") == 0) {
        return Hypothesis::TRIDIMENSIONAL;
      }
      raise("toString : unsupported modelling hypothesis");
    }  // end of fromString

    const char* toString(const Hypothesis h) {
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
    }  // end of toString

  }  // end of namespace behaviour

}  // end of namespace mgis
