/*!
 * \file   bindings/python/src/behaviour-module.cxx
 * \brief
 * \author Thomas Helfer
 * \date   31/10/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>

// forward declarations
void declareHypothesis(pybind11::module_&);
void declareVariable(pybind11::module_&);
void declareBehaviour(pybind11::module_&);
void declareBehaviourDescription(pybind11::module_&);
void declareState(pybind11::module_&);
void declareBehaviourData(pybind11::module_&);
void declareBehaviourDataView(pybind11::module_&);
void declareMaterialDataManager(pybind11::module_&);
void declareMaterialStateManager(pybind11::module_&);
void declareIntegrate(pybind11::module_&);
void declareFiniteStrainSupport(pybind11::module_&);

PYBIND11_MODULE(behaviour, m) {
  declareHypothesis(m);
  declareVariable(m);
  declareBehaviourDescription(m);
  declareBehaviour(m);
  declareState(m);
  declareBehaviourData(m);
  declareBehaviourDataView(m);
  declareMaterialDataManager(m);
  declareMaterialStateManager(m);
  declareIntegrate(m);
  declareFiniteStrainSupport(m);
}  // end of module behaviour
