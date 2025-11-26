/*!
 * \file   State.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <ostream>
#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Utilities/Markdown.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/StateView.hxx"
#include "MGIS/Behaviour/State.hxx"

namespace mgis::behaviour {

  State::State(State&&) = default;
  State::State(const State&) = default;

  State& State::operator=(State&& src) {
    if (&src.b != &this->b) {
      mgis::raise("State::operator=: unmatched behaviour");
    }
    if (&src != this) {
      this->gradients = std::move(src.gradients);
      this->thermodynamic_forces = std::move(src.thermodynamic_forces);
      this->material_properties = std::move(src.material_properties);
      this->internal_state_variables = std::move(src.internal_state_variables);
      this->stored_energy = src.stored_energy;
      this->dissipated_energy = src.dissipated_energy;
      this->external_state_variables = std::move(src.external_state_variables);
    }
    return *this;
  }  // end of State::operator=

  State& State::operator=(const State& src) {
    if (&src.b != &this->b) {
      mgis::raise("State::operator=: unmatched behaviour");
    }
    if (&src != this) {
      this->gradients = src.gradients;
      this->thermodynamic_forces = src.thermodynamic_forces;
      this->material_properties = src.material_properties;
      this->internal_state_variables = src.internal_state_variables;
      this->stored_energy = src.stored_energy;
      this->dissipated_energy = src.dissipated_energy;
      this->external_state_variables = src.external_state_variables;
    }
    return *this;
  }  // end of State::operator=

  State::State(const Behaviour& behaviour) : b(behaviour) {
    auto init = [this](std::vector<real>& values,
                       const std::vector<Variable> variables) {
      const auto s = getArraySize(variables, this->b.hypothesis);
      values.resize(s, real{0});
    };
    constexpr const auto zero = real{0};
    init(this->gradients, this->b.gradients);
    if ((this->b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (this->b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
      this->gradients[0] = this->gradients[1] = this->gradients[2] = real{1};
    }
    init(this->thermodynamic_forces, this->b.thermodynamic_forces);
    init(this->material_properties, this->b.mps);
    init(this->internal_state_variables, this->b.isvs);
    this->stored_energy = zero;
    this->dissipated_energy = zero;
    init(this->external_state_variables, this->b.esvs);
  }  // end of State::State

  void setGradient(State& s, const std::string_view n, const real v) {
    const auto& iv = getVariable(s.b.gradients, n);
    const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setGradient(s, o, v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setGradient(s, o, sz, v);
    }
  }  // end of setGradient

  void setGradient(State& s, const std::string_view n, const real* const v) {
    const auto& iv = getVariable(s.b.gradients, n);
    const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setGradient(s, o, *v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setGradient(s, o, sz, v);
    }
  }  // end of setGradient

  void setGradient(State& s, const size_type o, const real v) {
    s.gradients[o] = v;
  }  // end of setGradient

  void setGradient(State& s,
                   const size_type o,
                   const size_type n,
                   const real v) {
    auto p = s.gradients.begin() + o;
    std::fill(p, p + n, v);
  }  // end of setGradient

  void setGradient(State& s,
                   const size_type o,
                   const size_type n,
                   const real* const v) {
    auto p = s.gradients.begin() + o;
    std::copy(v, v + n, p);
  }  // end of setGradient

  real* getGradient(State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
    return getGradient(s, o);
  }  // end of getGradient

  const real* getGradient(const State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
    return getGradient(s, o);
  }  // end of getGradient

  real* getGradient(State& s, const size_type o) {
    return &(s.gradients[o]);
  }  // end of getGradient

  const real* getGradient(const State& s, const size_type o) {
    return &(s.gradients[o]);
  }  // end of getGradient

  void setThermodynamicForce(State& s, const std::string_view n, const real v) {
    const auto& iv = getVariable(s.b.thermodynamic_forces, n);
    const auto o =
        getVariableOffset(s.b.thermodynamic_forces, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setThermodynamicForce(s, o, v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setThermodynamicForce(s, o, sz, v);
    }
  }  // end of setThermodynamicForce

  void setThermodynamicForce(State& s,
                             const std::string_view n,
                             const real* const v) {
    const auto& iv = getVariable(s.b.thermodynamic_forces, n);
    const auto o =
        getVariableOffset(s.b.thermodynamic_forces, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setThermodynamicForce(s, o, *v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setThermodynamicForce(s, o, sz, v);
    }
  }  // end of setThermodynamicForce

  void setThermodynamicForce(State& s, const size_type o, const real v) {
    s.thermodynamic_forces[o] = v;
  }  // end of setThermodynamicForce

  void setThermodynamicForce(State& s,
                             const size_type o,
                             const size_type n,
                             const real v) {
    auto p = s.thermodynamic_forces.begin() + o;
    std::fill(p, p + n, v);
  }  // end of setThermodynamicForce

  void setThermodynamicForce(State& s,
                             const size_type o,
                             const size_type n,
                             const real* const v) {
    auto p = s.thermodynamic_forces.begin() + o;
    std::copy(v, v + n, p);
  }  // end of setThermodynamicForce

  real* getThermodynamicForce(State& s, const std::string_view n) {
    const auto o =
        getVariableOffset(s.b.thermodynamic_forces, n, s.b.hypothesis);
    return getThermodynamicForce(s, o);
  }  // end of getThermodynamicForce

  const real* getThermodynamicForce(const State& s, const std::string_view n) {
    const auto o =
        getVariableOffset(s.b.thermodynamic_forces, n, s.b.hypothesis);
    return getThermodynamicForce(s, o);
  }  // end of getThermodynamicForce

  real* getThermodynamicForce(State& s, const size_type o) {
    return &(s.thermodynamic_forces[o]);
  }  // end of getThermodynamicForce

  const real* getThermodynamicForce(const State& s, const size_type o) {
    return &(s.thermodynamic_forces[o]);
  }  // end of getThermodynamicForce

  void setMaterialProperty(State& s, const std::string_view n, const real v) {
    const auto o = getVariableOffset(s.b.mps, n, s.b.hypothesis);
    setMaterialProperty(s, o, v);
  }  // end of setMaterialProperty

  real* getMaterialProperty(State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.mps, n, s.b.hypothesis);
    return getMaterialProperty(s, o);
  }  // end of getMaterialProperty

  const real* getMaterialProperty(const State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.mps, n, s.b.hypothesis);
    return getMaterialProperty(s, o);
  }  // end of getMaterialProperty

  void setMaterialProperty(State& s, const size_type o, const real v) {
    s.material_properties[o] = v;
  }  // end of setMaterialProperty

  real* getMaterialProperty(State& s, const size_type o) {
    return &(s.material_properties[o]);
  }  // end of getMaterialProperty

  const real* getMaterialProperty(const State& s, const size_type o) {
    return &(s.material_properties[o]);
  }  // end of getMaterialProperty

  void setInternalStateVariable(State& s, const std::string_view n, const real v) {
    const auto& iv = getVariable(s.b.isvs, n);
    const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setInternalStateVariable(s, o, v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setInternalStateVariable(s, o, sz, v);
    }
  }  // end of setInternalStateVariable

  void setInternalStateVariable(State& s,
                                const std::string_view n,
                                const real* const v) {
    const auto& iv = getVariable(s.b.isvs, n);
    const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
    if (iv.type == Variable::SCALAR) {
      setInternalStateVariable(s, o, *v);
    } else {
      const auto sz = getVariableSize(iv, s.b.hypothesis);
      setInternalStateVariable(s, o, sz, v);
    }
  }  // end of setInternalStateVariable

  void setInternalStateVariable(State& s, const size_type o, const real v) {
    s.internal_state_variables[o] = v;
  }  // end of setInternalStateVariable

  void setInternalStateVariable(State& s,
                                const size_type o,
                                const size_type n,
                                const real v) {
    auto p = s.internal_state_variables.begin() + o;
    std::fill(p, p + n, v);
  }  // end of setInternalStateVariable

  void setInternalStateVariable(State& s,
                                const size_type o,
                                const size_type n,
                                const real* const v) {
    auto p = s.internal_state_variables.begin() + o;
    std::copy(v, v + n, p);
  }  // end of setInternalStateVariable

  real* getInternalStateVariable(State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
    return getInternalStateVariable(s, o);
  }  // end of getInternalStateVariable

  const real* getInternalStateVariable(const State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
    return getInternalStateVariable(s, o);
  }  // end of getInternalStateVariable

  real* getInternalStateVariable(State& s, const size_type o) {
    return &(s.internal_state_variables[o]);
  }  // end of getInternalStateVariable

  const real* getInternalStateVariable(const State& s, const size_type o) {
    return &(s.internal_state_variables[o]);
  }  // end of getInternalStateVariable

  void setExternalStateVariable(State& s, const std::string_view n, const real v) {
    const auto& ev = getVariable(s.b.esvs, n);
    if (ev.type != Variable::SCALAR) {
      mgis::raise("setExternalStateVariable: external state variable '" +
                  std::string{n} + "' is not a scalar");
    }
    const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
    setExternalStateVariable(s, o, v);
  }  // end of setExternalStateVariable

  void setExternalStateVariable(State& s,
                                const std::string_view n,
                                const mgis::span<const real> v) {
    const auto& ev = getVariable(s.b.esvs, n);
    const auto es = getVariableSize(ev, s.b.hypothesis);
    if (v.size() != es) {
      mgis::raise(
          "setExternalSateVariable: invalid number of values "
          "for external variable '" +
          std::string{n} +
          "' "
          "(" +
          std::to_string(v.size()) + " given, " + std::to_string(es) +
          "expected)");
    }
    const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
    setExternalStateVariable(s, o, v);
  }  // end of setExternalStateVariable

  void setExternalStateVariable(State& s, const size_type o, const real v) {
    s.external_state_variables[o] = v;
  }  // end of setExternalStateVariable

  void setExternalStateVariable(State& s,
                                const size_type o,
                                const mgis::span<const real> v) {
    std::copy(v.begin(), v.end(), s.external_state_variables.begin() + o);
  }  // end of setExternalStateVariable

  real* getExternalStateVariable(State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
    return getExternalStateVariable(s, o);
  }  // end of getExternalStateVariable

  const real* getExternalStateVariable(const State& s, const std::string_view n) {
    const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
    return getExternalStateVariable(s, o);
  }  // end of getExternalStateVariable

  real* getExternalStateVariable(State& s, const size_type o) {
    return &(s.external_state_variables[o]);
  }  // end of getExternalStateVariable

  const real* getExternalStateVariable(const State& s, const size_type o) {
    return &(s.external_state_variables[o]);
  }  // end of getExternalStateVariable

  StateView make_view(State& s) {
    auto get_ptr = [](std::vector<real>& v) -> real* {
      if (v.empty()) {
        return nullptr;
      }
      return &v[0];
    };  // end of get_ptr
    StateView v;
    v.mass_density = &(s.mass_density);
    v.gradients = get_ptr(s.gradients);
    v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
    v.material_properties = get_ptr(s.material_properties);
    v.internal_state_variables = get_ptr(s.internal_state_variables);
    v.stored_energy = &s.stored_energy;
    v.dissipated_energy = &s.dissipated_energy;
    v.external_state_variables = get_ptr(s.external_state_variables);
    return v;
  }  // end of make_view

  InitialStateView make_view(const State& s) {
    auto get_ptr = [](const std::vector<real>& v) -> const real* {
      if (v.empty()) {
        return nullptr;
      }
      return &v[0];
    };  // end of get_ptr
    InitialStateView v;
    v.mass_density = &(s.mass_density);
    v.gradients = get_ptr(s.gradients);
    v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
    v.material_properties = get_ptr(s.material_properties);
    v.internal_state_variables = get_ptr(s.internal_state_variables);
    v.stored_energy = &s.stored_energy;
    v.dissipated_energy = &s.dissipated_energy;
    v.external_state_variables = get_ptr(s.external_state_variables);
    return v;
  }  // end of make_view

  static void print_variables(std::ostream& os,
                              const Behaviour& b,
                              const std::vector<Variable>& variables,
                              const std::vector<mgis::real>& values) {
    auto o = mgis::size_type{};
    for (const auto& v : variables) {
      os << "- " << v.name << " (" << getVariableTypeAsString(v) << "): ";
      if (v.type == Variable::SCALAR) {
        if (values.size() < o) {
          mgis::raise("print_variables: invalid state initialisation");
        }
        os << values[o] << '\n';
        ++o;
      } else {
        const auto s = getVariableSize(v, b.hypothesis);
        if (values.size() < o + s) {
          mgis::raise("print_variables: invalid state initialisation");
        }
        os << '{';
        for (auto i = o; i != o + s;) {
          os << values[i];
          if (++i != o + s) {
            os << ", ";
          }
        }
        os << "}\n";
        o += s;
      }
    }
  }  // end of print_variables

  void print_markdown(std::ostream& os,
                      const Behaviour& b,
                      const State& s,
                      const mgis::size_type l) {
    if (!b.gradients.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1) << " Gradients\n\n";
      print_variables(os, b, b.gradients, s.gradients);
      os << '\n';
    }
    if (!b.thermodynamic_forces.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Thermodynamic forces\n\n";
      print_variables(os, b, b.thermodynamic_forces, s.thermodynamic_forces);
      os << '\n';
    }
    if (!b.mps.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Material properties\n\n";
      print_variables(os, b, b.mps, s.material_properties);
      os << '\n';
    }
    if (!b.isvs.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " Internal state variables\n\n";
      print_variables(os, b, b.isvs, s.internal_state_variables);
      os << '\n';
    }
    if (!b.esvs.empty()) {
      os << mgis::utilities::get_heading_signs(l + 1)
         << " External state variables\n\n";
      print_variables(os, b, b.esvs, s.external_state_variables);
      os << '\n';
    }
  }  // end of print_markdown

}  // end of namespace mgis::behaviour
