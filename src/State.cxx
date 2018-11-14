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

#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/StateView.hxx"
#include "MGIS/Behaviour/State.hxx"

namespace mgis {

  namespace behaviour {

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
        this->internal_state_variables =
            std::move(src.internal_state_variables);
        this->stored_energy = src.stored_energy;
        this->dissipated_energy = src.dissipated_energy;
        this->external_state_variables =
            std::move(src.external_state_variables);
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
      init(this->thermodynamic_forces, this->b.thermodynamic_forces);
      init(this->material_properties, this->b.mps);
      init(this->internal_state_variables, this->b.isvs);
      this->stored_energy = zero;
      this->dissipated_energy = zero;
      init(this->external_state_variables, this->b.esvs);
    }  // end of State::State

    void setGradient(State& s, const string_view n, const real v) {
      const auto& iv = getVariable(s.b.gradients, n);
      const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setGradient(s, o, v);
      } else {
        const auto sz = getVariableSize(iv, s.b.hypothesis);
        setGradient(s, o, sz, v);
      }
    }  // end of setGradient

    void setGradient(State& s, const string_view n, const real* const v) {
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

    real* getGradient(State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
      return getGradient(s, o);
    }  // end of getGradient

    const real* getGradient(const State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.gradients, n, s.b.hypothesis);
      return getGradient(s, o);
    }  // end of getGradient

    real* getGradient(State& s, const size_type o) {
      return &(s.gradients[o]);
    }  // end of getGradient

    const real* getGradient(const State& s, const size_type o) {
      return &(s.gradients[o]);
    }  // end of getGradient

    void setThermodynamicForce(State& s, const string_view n, const real v) {
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
                               const string_view n,
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

    real* getThermodynamicForce(State& s, const string_view n) {
      const auto o =
          getVariableOffset(s.b.thermodynamic_forces, n, s.b.hypothesis);
      return getThermodynamicForce(s, o);
    }  // end of getThermodynamicForce

    const real* getThermodynamicForce(const State& s, const string_view n) {
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

    void setMaterialProperty(State& s, const string_view n, const real v) {
      const auto o = getVariableOffset(s.b.mps, n, s.b.hypothesis);
      setMaterialProperty(s, o, v);
    }  // end of setMaterialProperty

    real* getMaterialProperty(State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.mps, n, s.b.hypothesis);
      return getMaterialProperty(s, o);
    }  // end of getMaterialProperty

    const real* getMaterialProperty(const State& s, const string_view n) {
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

    void setInternalStateVariable(State& s, const string_view n, const real v) {
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
                                  const string_view n,
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

    real* getInternalStateVariable(State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
      return getInternalStateVariable(s, o);
    }  // end of getInternalStateVariable

    const real* getInternalStateVariable(const State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.isvs, n, s.b.hypothesis);
      return getInternalStateVariable(s, o);
    }  // end of getInternalStateVariable

    real* getInternalStateVariable(State& s, const size_type o) {
      return &(s.internal_state_variables[o]);
    }  // end of getInternalStateVariable

    const real* getInternalStateVariable(const State& s, const size_type o) {
      return &(s.internal_state_variables[o]);
    }  // end of getInternalStateVariable

    void setExternalStateVariable(State& s, const string_view n, const real v) {
      const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
      setExternalStateVariable(s, o, v);
    }  // end of setExternalStateVariable

    real* getExternalStateVariable(State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
      return getExternalStateVariable(s, o);
    }  // end of getExternalStateVariable

    const real* getExternalStateVariable(const State& s, const string_view n) {
      const auto o = getVariableOffset(s.b.esvs, n, s.b.hypothesis);
      return getExternalStateVariable(s, o);
    }  // end of getExternalStateVariable

    void setExternalStateVariable(State& s, const size_type o, const real v) {
      s.external_state_variables[o] = v;
    }  // end of setExternalStateVariable

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
      v.gradients = get_ptr(s.gradients);
      v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
      v.material_properties = get_ptr(s.material_properties);
      v.internal_state_variables = get_ptr(s.internal_state_variables);
      v.stored_energy = &s.stored_energy;
      v.dissipated_energy = &s.dissipated_energy;
      v.external_state_variables = get_ptr(s.external_state_variables);
      return v;
    }  // end of make_view

  }  // end of namespace behaviour

}  // end of namespace mgis
