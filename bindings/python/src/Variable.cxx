/*!
 * \file   Variable.cxx
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
#include <pybind11/stl.h>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"

static const char* Variable_getType(const mgis::behaviour::Variable& v) {
  using mgis::behaviour::Variable;
  switch (v.type) {
    case Variable::SCALAR:
      return "Scalar";
    case Variable::VECTOR:
      return "Vector";
    case Variable::STENSOR:
      return "Stensor";
    case Variable::TENSOR:
      return "Tensor";
    default:
      mgis::raise("Variable_getType: unsupported type");
  }
  return "";
}  // end of Variable_getType

// mgis::string_view is not exposed
static mgis::size_type getVariableOffsetByString(
    const std::vector<mgis::behaviour::Variable>& vs,
    const std::string& n,
    const mgis::behaviour::Hypothesis h) {
  return mgis::behaviour::getVariableOffset(vs, n, h);
}  // getVariableOffsetByString

// mgis::string_view is not exposed
static const mgis::behaviour::Variable& getVariableByString(
    const std::vector<mgis::behaviour::Variable>& vs, const std::string& n) {
  return mgis::behaviour::getVariable(vs, n);
}  // getVariableByString

// mgis::string_view is not exposed
static mgis::size_type getVariableSizeByString(
    const std::vector<mgis::behaviour::Variable>& vs,
    const std::string& n,
    const mgis::behaviour::Hypothesis h) {
  return mgis::behaviour::getVariableSize(mgis::behaviour::getVariable(vs, n),
                                          h);
}  // getVariableSizeByString

// forward declaration
void declareVariable(pybind11::module_&);

void declareVariable(pybind11::module_& m) {
  using mgis::behaviour::Variable;
  // wrapping the Variable::Type enum
  pybind11::enum_<Variable::Type>(m, "VariableType")
      .value("Scalar", Variable::SCALAR)
      .value("SCALAR", Variable::SCALAR)
      .value("Vector", Variable::VECTOR)
      .value("Vector1D", Variable::VECTOR_1D)
      .value("Vector2D", Variable::VECTOR_2D)
      .value("Vector3D", Variable::VECTOR_3D)
      .value("VECTOR", Variable::VECTOR)
      .value("VECTOR1D", Variable::VECTOR_1D)
      .value("VECTOR2D", Variable::VECTOR_2D)
      .value("VECTOR3D", Variable::VECTOR_3D)
      .value("Stensor", Variable::STENSOR)
      .value("Stensor1D", Variable::STENSOR_1D)
      .value("Stensor2D", Variable::STENSOR_2D)
      .value("Stensor3D", Variable::STENSOR_3D)
      .value("STENSOR", Variable::STENSOR)
      .value("STENSOR1D", Variable::STENSOR_1D)
      .value("STENSOR2D", Variable::STENSOR_2D)
      .value("STENSOR3D", Variable::STENSOR_3D)
      .value("Tensor", Variable::TENSOR)
      .value("Tensor1D", Variable::TENSOR_1D)
      .value("Tensor2D", Variable::TENSOR_2D)
      .value("Tensor3D", Variable::TENSOR_3D)
      .value("TENSOR", Variable::TENSOR)
      .value("TENSOR1D", Variable::TENSOR_1D)
      .value("TENSOR2D", Variable::TENSOR_2D)
      .value("TENSOR3D", Variable::TENSOR_3D)
      .value("HIGHER_ORDER_TENSOR", Variable::HIGHER_ORDER_TENSOR)
      .value("HigherOrderTensor", Variable::HIGHER_ORDER_TENSOR)
      .value("ARRAY", Variable::ARRAY)
      .value("Array", Variable::ARRAY);
  // wrapping the Variable class
  pybind11::class_<Variable>(m, "Variable")
      .def_readonly("name", &Variable::name, "the name of the variable")
      .def_readonly("type", &Variable::type, "the type of the variable.")
      .def("getType", Variable_getType,
           "the type of the variable. "
           "Possible values are `Scalar`, `Vector`, `Stensor`, `Tensor`");
  // free functions
  m.def("getVariable", getVariableByString,
        pybind11::return_value_policy::reference);
  m.def("getVariableSize", &mgis::behaviour::getVariableSize);
  m.def("getVariableSize", getVariableSizeByString);
  m.def("getArraySize", &mgis::behaviour::getArraySize);
  m.def("getVariableOffset", getVariableOffsetByString);
  m.def("getVariableTypeSymbolicRepresentation",
        &mgis::behaviour::getVariableTypeSymbolicRepresentation);

}  // end of declareVariable
