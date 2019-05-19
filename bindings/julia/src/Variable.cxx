/*!
 * \file   bindings/julia/src/Variable.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <stdexcept>
#include <jlcxx/jlcxx.hpp>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareVariable(jlcxx::Module& m) {
  using mgis::behaviour::Variable;
  m.add_bits<Variable::Type>("VariableType");
  m.set_const("Scalar", Variable::SCALAR);
  m.set_const("SCALAR", Variable::SCALAR);
  m.set_const("Vector", Variable::VECTOR);
  m.set_const("VECTOR", Variable::VECTOR);
  m.set_const("Stensor", Variable::STENSOR);
  m.set_const("STENSOR", Variable::STENSOR);
  m.set_const("Tensor", Variable::TENSOR);
  m.set_const("TENSOR", Variable::TENSOR);
  m.add_type<Variable>("Variable")
      .method("get_name", [](const Variable& v) { return v.name; })
      .method("get_type", [](const Variable& v) { return v.type; });
  mgis::julia::expose_std_vector<Variable>(m, "VariablesVector");
}  // end of declareVariable

