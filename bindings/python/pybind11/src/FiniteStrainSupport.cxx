/*!
 * \file   FiniteStrainSupport.cxx
 * \brief
 * \author Thomas Helfer
 * \date   25/01/2019
 * \copyright (C) Copyright Thomas Helfer 2018-2019.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/FiniteStrainSupport.hxx"

// forward declaration
void declareFiniteStrainSupport(pybind11::module_&);

static void py_convertFiniteStrainStress(
    pybind11::object o,
    const mgis::behaviour::MaterialDataManager& m,
    const mgis::behaviour::FiniteStrainStress t) {
  auto s = mgis::python::mgis_convert_to_span(o);
  mgis::behaviour::convertFiniteStrainStress(s, m, t);
}  // end of py_py_convertFiniteStrainStress

static void py_convertFiniteStrainTangentOperator(
    pybind11::object o,
    const mgis::behaviour::MaterialDataManager& m,
    const mgis::behaviour::FiniteStrainTangentOperator t) {
  auto K = mgis::python::mgis_convert_to_span(o);
  mgis::behaviour::convertFiniteStrainTangentOperator(K, m, t);
}  // end of py_py_convertFiniteStrainTangentOperator

void declareFiniteStrainSupport(pybind11::module_& m) {
  using namespace mgis::behaviour;
  pybind11::enum_<FiniteStrainStress>(m, "FiniteStrainStress")
      .value("PK1", FiniteStrainStress::PK1)
      .value("FirstPiolaKirchhoffStress", FiniteStrainStress::PK1);
  pybind11::enum_<FiniteStrainTangentOperator>(m, "FiniteStrainTangentOperator")
      .value("DPK1_DF", FiniteStrainTangentOperator::DPK1_DF);

  m.def("convertFiniteStrainStress", py_convertFiniteStrainStress);
  m.def("convertFiniteStrainTangentOperator",
        py_convertFiniteStrainTangentOperator);
}  // end of declareFiniteStrainSupport
