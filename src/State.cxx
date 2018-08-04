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
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/State.hxx"

namespace mgis {

  namespace behaviour {

    State::State() = default;
    State::State(State&&) = default;
    State::State(const State&) = default;
    State& State::operator=(State&&) = default;
    State& State::operator=(const State&) = default;

    State::State(const Behaviour& b){
      auto init = [&b](std::vector<real>& values,
                       const std::vector<Variable> variables) {
        constexpr const auto zero = real{0};
        const auto s = getArraySize(variables, b.hypothesis);
        values.resize(s, zero);
      };
      init(this->gradients, b.gradients);
      init(this->thermodynamic_forces, b.thermodynamic_forces);
      init(this->material_properties, b.mps);
      init(this->internal_state_variables, b.isvs);
      init(this->external_state_variables, b.esvs);
    }  // end of State::State

    void setGradient(State& s,
                     const Behaviour& b,
                     const string_view n,
                     const real v) {
      const auto& iv = getVariable(b.gradients, n);
      const auto o = getVariableOffset(b.gradients, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setGradient(s, o, v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
        setGradient(s, o, sz, v);
      }
    }  // end of setGradient

    void setGradient(State& s,
                     const Behaviour& b,
                     const string_view n,
                     const real* const v) {
      const auto& iv = getVariable(b.gradients, n);
      const auto o = getVariableOffset(b.gradients, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setGradient(s, o, *v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
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

    real* getGradient(State& s, const Behaviour& b, const string_view n) {
      const auto o = getVariableOffset(b.gradients, n, b.hypothesis);
      return getGradient(s, o);
    } // end of getGradient

    const real* getGradient(const State& s,
                            const Behaviour& b,
                            const string_view n) {
      const auto o = getVariableOffset(b.gradients, n, b.hypothesis);
      return getGradient(s, o);
    } // end of getGradient

    real* getGradient(State& s, const size_type o) {
      return &(s.gradients[o]);
    } // end of getGradient

    const real* getGradient(const State& s, const size_type o) {
      return &(s.gradients[o]);
    } // end of getGradient

    void setThermodynamicForce(State& s,
                               const Behaviour& b,
                               const string_view n,
                               const real v) {
      const auto& iv = getVariable(b.thermodynamic_forces, n);
      const auto o = getVariableOffset(b.thermodynamic_forces, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setThermodynamicForce(s, o, v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
        setThermodynamicForce(s, o, sz, v);
      }
    }  // end of setThermodynamicForce

    void setThermodynamicForce(State& s,
                               const Behaviour& b,
                               const string_view n,
                               const real* const v) {
      const auto& iv = getVariable(b.thermodynamic_forces, n);
      const auto o = getVariableOffset(b.thermodynamic_forces, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setThermodynamicForce(s, o, *v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
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

    real* getThermodynamicForce(State& s,
                                const Behaviour& b,
                                const string_view n) {
      const auto o = getVariableOffset(b.thermodynamic_forces, n, b.hypothesis);
      return getThermodynamicForce(s, o);
    } // end of getThermodynamicForce

    const real* getThermodynamicForce(const State& s,
                                      const Behaviour& b,
                                      const string_view n) {
      const auto o = getVariableOffset(b.thermodynamic_forces, n, b.hypothesis);
      return getThermodynamicForce(s, o);
    } // end of getThermodynamicForce

    real* getThermodynamicForce(State& s, const size_type o) {
      return &(s.thermodynamic_forces[o]);
    } // end of getThermodynamicForce

    const real* getThermodynamicForce(const State& s, const size_type o) {
      return &(s.thermodynamic_forces[o]);
    } // end of getThermodynamicForce

    void setMaterialProperty(State& s,
                             const Behaviour& b,
                             const string_view n,
                             const real v) {
      const auto o = getVariableOffset(b.mps, n, b.hypothesis);
      setMaterialProperty(s, o, v);
    }  // end of setMaterialProperty

    real* getMaterialProperty(State& s,
                              const Behaviour& b,
                              const string_view n) {
      const auto o = getVariableOffset(b.mps, n, b.hypothesis);
      return getMaterialProperty(s, o);
    }  // end of getMaterialProperty

    const real* getMaterialProperty(const State& s,
                                    const Behaviour& b,
                                    const string_view n) {
      const auto o = getVariableOffset(b.mps, n, b.hypothesis);
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

    void setInternalStateVariable(State& s,
                                  const Behaviour& b,
                                  const string_view n,
                                  const real v) {
      const auto& iv = getVariable(b.isvs, n);
      const auto o = getVariableOffset(b.isvs, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setInternalStateVariable(s, o, v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
        setInternalStateVariable(s, o, sz, v);
      }
    }  // end of setInternalStateVariable

    void setInternalStateVariable(State& s,
                                  const Behaviour& b,
                                  const string_view n,
                                  const real* const v) {
      const auto& iv = getVariable(b.isvs, n);
      const auto o = getVariableOffset(b.isvs, n, b.hypothesis);
      if (iv.type == Variable::SCALAR) {
        setInternalStateVariable(s, o, *v);
      } else {
        const auto sz = getVariableSize(iv, b.hypothesis);
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

    real* getInternalStateVariable(State& s,
                                   const Behaviour& b,
                                   const string_view n) {
      const auto o = getVariableOffset(b.isvs, n, b.hypothesis);
      return getInternalStateVariable(s, o);
    } // end of getInternalStateVariable

    const real* getInternalStateVariable(const State& s,
                                         const Behaviour& b,
                                         const string_view n) {
      const auto o = getVariableOffset(b.isvs, n, b.hypothesis);
      return getInternalStateVariable(s, o);
    } // end of getInternalStateVariable

    real* getInternalStateVariable(State& s, const size_type o) {
      return &(s.internal_state_variables[o]);
    } // end of getInternalStateVariable

    const real* getInternalStateVariable(const State& s, const size_type o) {
      return &(s.internal_state_variables[o]);
    } // end of getInternalStateVariable

    void setExternalStateVariable(State& s,
                                  const Behaviour& b,
                                  const string_view n,
                                  const real v) {
      const auto o = getVariableOffset(b.esvs, n, b.hypothesis);
      setExternalStateVariable(s, o, v);
    }  // end of setExternalStateVariable

    real* getExternalStateVariable(State& s,
                                   const Behaviour& b,
                                   const string_view n) {
      const auto o = getVariableOffset(b.esvs, n, b.hypothesis);
      return getExternalStateVariable(s, o);
    }  // end of getExternalStateVariable

    const real* getExternalStateVariable(const State& s,
                                         const Behaviour& b,
                                         const string_view n) {
      const auto o = getVariableOffset(b.esvs, n, b.hypothesis);
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

  }  // end of namespace behaviour

}  // end of namespace mgis
