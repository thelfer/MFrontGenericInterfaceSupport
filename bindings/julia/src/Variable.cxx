/*!
 * \file   bindings/julia/src/Variable.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   14/05/2019
 */

#include <cstdint>
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
      .method("get_type", [](const Variable& v) { return v.type; })
      .method("get_variable_size",
              [](const Variable& v, const mgis::behaviour::Hypothesis h) {
                return mgis::behaviour::getVariableSize(v, h);
              })
      .method("get_size",
              [](const Variable& v, const mgis::behaviour::Hypothesis h) {
                return mgis::behaviour::getVariableSize(v, h);
              });

  mgis::julia::expose_std_vector<Variable>(m, "VariablesVector");
  m.method("contains",
           [](std::vector<Variable>& v, const std::string& n) {
             return mgis::behaviour::contains(v, n);
           });
  m.method("get_variable",
           [](const std::vector<Variable>& v, const std::string& n) {
             return mgis::behaviour::getVariable(v, n);
           });
  m.method("get_array_size", [](const std::vector<Variable>& v,
                                const mgis::behaviour::Hypothesis h) {
    return mgis::behaviour::getArraySize(v, h);
  });
  m.method("get_variable_offset",
           [](const std::vector<Variable>& v, const std::string& n,
              const mgis::behaviour::Hypothesis h) -> std::int64_t {
             return mgis::behaviour::getVariableOffset(v, n, h) + 1;
           });

}  // end of declareVariable

