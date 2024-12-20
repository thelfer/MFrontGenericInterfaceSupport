/*!
 * \file   Hypothesis.cxx
 * \brief
 * \author Thomas Helfer
 * \date   06/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>
#include "MGIS/Behaviour/Hypothesis.hxx"

// forward declaration
void declareHypothesis(pybind11::module_&);

void declareHypothesis(pybind11::module_& m) {
  using mgis::behaviour::Hypothesis;
  pybind11::enum_<mgis::behaviour::Hypothesis>(m, "Hypothesis")
      .value("AXISYMMETRICALGENERALISEDPLANESTRAIN",
             Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN)
      .value("AxisymmetricalGeneralisedPlaneStrain",
             Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN)
      .value("AxisymmetricalGeneralisedPlaneStress",
             Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)
      .value("AXISYMMETRICALGENERALISEDPLANESTRESS",
             Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)
      .value("Axisymmetrical", Hypothesis::AXISYMMETRICAL)
      .value("AXISYMMETRICAL", Hypothesis::AXISYMMETRICAL)
      .value("PlaneStress", Hypothesis::PLANESTRESS)
      .value("PLANESTRESS", Hypothesis::PLANESTRESS)
      .value("PlaneStrain", Hypothesis::PLANESTRAIN)
      .value("PLANESTRAIN", Hypothesis::PLANESTRAIN)
      .value("GeneralisedPlaneStrain", Hypothesis::GENERALISEDPLANESTRAIN)
      .value("GENERALISEDPLANESTRAIN", Hypothesis::GENERALISEDPLANESTRAIN)
      .value("Tridimensional", Hypothesis::TRIDIMENSIONAL)
      .value("TRIDIMENSIONAL", Hypothesis::TRIDIMENSIONAL);
}  // end of declareHypothesis
