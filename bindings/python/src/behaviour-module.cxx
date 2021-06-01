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

#include <boost/python/module.hpp>
#include "MGIS/Python/NumPySupport.hxx"

// forward declarations
void declareHypothesis();
void declareVariable();
void declareBehaviour();
void declareState();
void declareBehaviourData();
void declareBehaviourDataView();
void declareMaterialDataManager();
void declareMaterialStateManager();
void declareIntegrate();
void declareFiniteStrainSupport();

BOOST_PYTHON_MODULE(behaviour) {
  mgis::python::initializeNumPy();
  declareHypothesis();
  declareVariable();
  declareBehaviour();
  declareState();
  declareBehaviourData();
  declareBehaviourDataView();
  declareMaterialDataManager();
  declareMaterialStateManager();
  declareIntegrate();
  declareFiniteStrainSupport();
}  // end of module behaviour
