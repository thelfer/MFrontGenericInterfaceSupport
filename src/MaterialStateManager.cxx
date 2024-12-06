/*!
 * \file   MaterialStateManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/08/2018
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
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis::behaviour {

  MaterialStateManager::MaterialStateManager(const Behaviour& behaviour,
                                             const size_type s)
      : gradients_stride(
            getArraySize(behaviour.gradients, behaviour.hypothesis)),
        thermodynamic_forces_stride(
            getArraySize(behaviour.thermodynamic_forces, behaviour.hypothesis)),
        internal_state_variables_stride(
            getArraySize(behaviour.isvs, behaviour.hypothesis)),
        n(s),
        b(behaviour) {
    auto init = [this](mgis::span<mgis::real>& view,
                       std::vector<mgis::real>& values, const size_type vs) {
      constexpr const auto zero = real{0};
      values.resize(this->n * vs, zero);
      view = mgis::span<mgis::real>(values);
    };
    init(this->gradients, this->gradients_values, this->gradients_stride);
    if ((this->b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (this->b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY)) {
      for (size_type i = 0; i != this->n; ++i) {
        auto F = this->gradients.subspan(i * gradients_stride,  //
                                         gradients_stride);
        F[0] = F[1] = F[2] = real{1};
      }
    }
    init(this->thermodynamic_forces, this->thermodynamic_forces_values,
         this->thermodynamic_forces_stride);
    init(this->internal_state_variables, this->internal_state_variables_values,
         this->internal_state_variables_stride);
    if (b.computesStoredEnergy) {
      init(this->stored_energies, this->stored_energies_values, 1u);
    }
    if (b.computesDissipatedEnergy) {
      init(this->dissipated_energies, this->dissipated_energies_values, 1u);
    }
  }  // end of MaterialStateManager::MaterialStateManager

  MaterialStateManager::MaterialStateManager(
      const Behaviour& behaviour,
      const size_type s,
      const MaterialStateManagerInitializer& i)
      : gradients_stride(
            getArraySize(behaviour.gradients, behaviour.hypothesis)),
        thermodynamic_forces_stride(
            getArraySize(behaviour.thermodynamic_forces, behaviour.hypothesis)),
        internal_state_variables_stride(
            getArraySize(behaviour.isvs, behaviour.hypothesis)),
        n(s),
        b(behaviour) {
    auto init = [this](mgis::span<mgis::real>& view,
                       std::vector<mgis::real>& values,
                       const mgis::span<mgis::real>& evalues,
                       const size_type vs, const char* const vn) {
      if (evalues.empty()) {
        constexpr const auto zero = real{0};
        values.resize(this->n * vs, zero);
        view = mgis::span<mgis::real>(values);
      } else {
        if (static_cast<size_type>(evalues.size()) != this->n * vs) {
          mgis::raise(
              "MaterialStateManager::MaterialStateManager: "
              "the memory associated with the " +
              std::string(vn) + " has not been allocated properly");
        }
        view = evalues;
      }
    };
    init(this->gradients, this->gradients_values, i.gradients,
         this->gradients_stride, "gradients");
    if ((this->b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (this->b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY) &&
        (i.gradients.empty())) {
      for (size_type ig = 0; ig != this->n; ++ig) {
        auto F = this->gradients.subspan(ig * gradients_stride,  //
                                         gradients_stride);
        F[0] = F[1] = F[2] = real{1};
      }
    }
    init(this->thermodynamic_forces, this->thermodynamic_forces_values,
         i.thermodynamic_forces, this->thermodynamic_forces_stride,
         "thermodynamic forces");
    init(this->internal_state_variables, this->internal_state_variables_values,
         i.internal_state_variables, this->internal_state_variables_stride,
         "internal state variables");
    if (b.computesStoredEnergy) {
      init(this->stored_energies, this->stored_energies_values,
           i.stored_energies, 1u, "stored energies");
    } else {
      if (!i.stored_energies.empty()) {
        mgis::raise(
            "MaterialStateManager::MaterialStateManager: "
            "stored energies shall not have been allocated as the behaviour "
            "don't compute the stored energy");
      }
    }
    if (b.computesDissipatedEnergy) {
      init(this->dissipated_energies, this->dissipated_energies_values,
           i.dissipated_energies, 1u, "dissipated energies");
    } else {
      if (!i.dissipated_energies.empty()) {
        mgis::raise(
            "MaterialStateManager::MaterialStateManager: "
            "dissipated energies shall not have been allocated as the "
            "behaviour don't compute the dissipated energy");
      }
    }
  }  // end of MaterialStateManager::MaterialStateManager

  MaterialStateManager::~MaterialStateManager() = default;

  static MaterialStateManager::FieldHolder& getFieldHolder(
      std::map<std::string, MaterialStateManager::FieldHolder>& m,
      const std::string_view& n) {
    // #if __cplusplus > 201103L
    //       return m[n];
    // #else  /* __cplusplus > 201103L */
    return m[std::string{n}];
    // #endif /* __cplusplus > 201103L */
  }  // end of getFieldHolder

  //   static std::map<std::string, MaterialStateManager::FieldHolder>::iterator
  //   getFieldHolderIterator(
  //       std::map<std::string, MaterialStateManager::FieldHolder>& m,
  //       const std::string_view& n) {
  //     // #if __cplusplus > 201103L
  //     //       return m.find(n);
  //     // #else  /* __cplusplus > 201103L */
  //     return m.find(n.to_string());
  //     // #endif /* __cplusplus > 201103L */
  //   }  // end of getFieldHolder

  static std::map<std::string,
                  MaterialStateManager::FieldHolder>::const_iterator
  getFieldHolderIterator(
      const std::map<std::string, MaterialStateManager::FieldHolder>& m,
      const std::string_view& n) {
    // #if __cplusplus > 201103L
    //       return m.find(n);
    // #else  /* __cplusplus > 201103L */
    return m.find(std::string{n});
    // #endif /* __cplusplus > 201103L */
  }  // end of getFieldHolder

  void setMaterialProperty(MaterialStateManager& m,
                           const std::string_view& n,
                           const real v) {
    const auto mp = getVariable(m.b.mps, n);
    mgis::raise_if(mp.type != Variable::SCALAR,
                   "setMaterialProperty: "
                   "invalid material property "
                   "(only scalar material property is supported)");
    getFieldHolder(m.material_properties, n) = v;
  }  // end of setMaterialProperty

  MGIS_EXPORT void setMaterialProperty(
      MaterialStateManager& m,
      const std::string_view& n,
      const mgis::span<real>& v,
      const MaterialStateManager::StorageMode s) {
    const auto mp = getVariable(m.b.mps, n);
    mgis::raise_if(mp.type != Variable::SCALAR,
                   "setMaterialProperty: "
                   "invalid material property "
                   "(only scalar material property is supported)");
    mgis::raise_if(static_cast<mgis::size_type>(v.size()) != m.n,
                   "setMaterialProperty: invalid number of values "
                   "(does not match the number of integration points)");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      getFieldHolder(m.material_properties, n) =
          std::vector<real>{v.begin(), v.end()};
    } else {
      getFieldHolder(m.material_properties, n) = v;
    }
  }  // end of setMaterialProperty

  bool isMaterialPropertyDefined(const MaterialStateManager& m,
                                 const std::string_view& n) {
    const auto p = getFieldHolderIterator(m.material_properties, n);
    return p != m.material_properties.end();
  }  // end of isMaterialPropertyDefined

  bool isMaterialPropertyUniform(const MaterialStateManager& m,
                                 const std::string_view& n) {
    const auto p = getFieldHolderIterator(m.material_properties, n);
    if (p == m.material_properties.end()) {
      mgis::raise(
          "isMaterialPropertyUniform: "
          "no material property named '" +
          std::string{n} + "' defined");
    }
    return std::holds_alternative<real>(p->second);
  }  // end of isMaterialPropertyUniform

  void setMassDensity(MaterialStateManager& m, const real v) {
    m.mass_density = v;
  }  // end of setMassDensity

  MGIS_EXPORT void setMassDensity(
      MaterialStateManager& m,
      const mgis::span<real>& v,
      const MaterialStateManager::StorageMode s) {
    mgis::raise_if(static_cast<mgis::size_type>(v.size()) != m.n,
                   "setMassDensity: invalid number of values "
                   "(does not match the number of integration points)");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      m.mass_density = std::vector<real>{v.begin(), v.end()};
    } else {
      m.mass_density = v;
    }
  }  // end of setMassDensity

  bool isMassDensityDefined(const MaterialStateManager& m) {
    return m.mass_density.has_value();
  }  // end of isMassDensityDefined

  bool isMassDensityUniform(const MaterialStateManager& m) {
    if (!isMassDensityDefined(m)) {
      mgis::raise("isMassDensityUniform: the mass density is undefined");
    }
    return std::holds_alternative<real>(*(m.mass_density));
  }  // end of isMassDensityUniform

  void setExternalStateVariable(MaterialStateManager& m,
                                const std::string_view& n,
                                const real v) {
    const auto esv = getVariable(m.b.esvs, n);
    mgis::raise_if(esv.type != Variable::SCALAR,
                   "setExternalStateVariable: "
                   "invalid external state variable "
                   "(only scalar external state variable is supported)");
    getFieldHolder(m.external_state_variables, n) = v;
  }  // end of setExternalStateVariable

  MGIS_EXPORT void setExternalStateVariable(
      MaterialStateManager& m,
      const std::string_view& n,
      const mgis::span<real>& v,
      const MaterialStateManager::StorageMode s) {
    const auto esv = getVariable(m.b.esvs, n);
    const auto vs = getVariableSize(esv, m.b.hypothesis);
    mgis::raise_if(((static_cast<mgis::size_type>(v.size()) != m.n * vs) &&
                    (static_cast<mgis::size_type>(v.size()) != vs)),
                   "setExternalStateVariable: invalid number of values");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      if (v.size() == 1u) {
        getFieldHolder(m.external_state_variables, n) = v[0];
      } else {
        getFieldHolder(m.external_state_variables, n) =
            std::vector<real>{v.begin(), v.end()};
      }
    } else {
      getFieldHolder(m.external_state_variables, n) = v;
    }
  }  // end of setExternalStateVariable

  bool isExternalStateVariableDefined(const MaterialStateManager& m,
                                      const std::string_view& n) {
    const auto p = getFieldHolderIterator(m.external_state_variables, n);
    return p != m.external_state_variables.end();
  }  // end of isExternalStateVariableDefined

  bool isExternalStateVariableUniform(const MaterialStateManager& m,
                                      const std::string_view& n) {
    const auto p = getFieldHolderIterator(m.external_state_variables, n);
    mgis::raise_if(p == m.external_state_variables.end(),
                   "isExternalStateVariableUniform: "
                   "no external state variable named '" +
                       std::string{n} + "' defined");
    return std::holds_alternative<real>(p->second);
  }  // end of isExternalStateVariableUniform

  void updateValues(MaterialStateManager& o, const MaterialStateManager& i) {
    auto check_size = [](const mgis::size_type s1, const mgis::size_type s2) {
      if (s1 != s2) {
        mgis::raise(
            "mgis::behaviour::updateValues: "
            "arrays' size does not match");
      }
    };  // end of check_size
    auto update_span = [&check_size](mgis::span<real>& to,
                                     const mgis::span<const real>& from) {
      check_size(from.size(), to.size());
      std::copy(from.begin(), from.end(), to.begin());
    };  // end update_span
    auto update_field_holder =
        [&check_size](MaterialStateManager::FieldHolder& to,
                      const MaterialStateManager::FieldHolder& from) {
          if (std::holds_alternative<mgis::real>(from)) {
            to = std::get<mgis::real>(from);
          } else if (std::holds_alternative<std::vector<mgis::real>>(from)) {
            const auto& from_v = std::get<std::vector<mgis::real>>(from);
            if (std::holds_alternative<mgis::span<mgis::real>>(to)) {
              // reuse existing memory
              auto& to_v = std::get<mgis::span<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else if (std::holds_alternative<std::vector<mgis::real>>(to)) {
              // reuse existing memory
              auto& to_v = std::get<std::vector<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else {
              // to contains a real value, so overwrite it with a new vector
              to = std::get<std::vector<mgis::real>>(from);
            }
          } else {
            const auto from_v = std::get<mgis::span<mgis::real>>(from);
            if (std::holds_alternative<mgis::span<mgis::real>>(to)) {
              // reuse existing memory
              auto to_v = std::get<mgis::span<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else if (std::holds_alternative<std::vector<mgis::real>>(to)) {
              // reuse existing memory
              auto to_v = std::get<std::vector<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else {
              to = from_v;
            }
          }
        };  // end of update_field_holder
    auto check_mps =
        [](const Behaviour& b,
           const std::map<std::string, MaterialStateManager::FieldHolder>&
               mps) {
          for (const auto& mp : mps) {
            auto find_mp = [&mp](const Variable& d) {
              return mp.first == d.name;
            };
            if (std::find_if(b.mps.begin(), b.mps.end(), find_mp) ==
                b.mps.end()) {
              mgis::raise(
                  "mgis::behaviour::updateValues: "
                  "material property '" +
                  mp.first +
                  "' defined in the material state manager is not defined "
                  " by the behaviour");
            }
          }
        };
    if (&i.b != &o.b) {
      mgis::raise(
          "mgis::behaviour::updateValues: the material state managers "
          "do not holds the same behaviour");
    }
    check_mps(o.b, i.material_properties);
    check_mps(o.b, o.material_properties);
    update_span(o.gradients, i.gradients);
    update_span(o.thermodynamic_forces, i.thermodynamic_forces);
    update_span(o.internal_state_variables, i.internal_state_variables);
    update_span(o.stored_energies, i.stored_energies);
    update_span(o.dissipated_energies, i.dissipated_energies);
    auto pmp = o.material_properties.begin();
    while (pmp != o.material_properties.end()) {
      if (i.material_properties.count(pmp->first) == 0) {
        pmp = o.material_properties.erase(pmp);
      } else {
        ++pmp;
      }
    }
    for (const auto& mp : i.material_properties) {
      update_field_holder(o.material_properties[mp.first], mp.second);
    }
    auto pev = o.external_state_variables.begin();
    while (pev != o.external_state_variables.end()) {
      if (i.external_state_variables.count(pev->first) == 0) {
        pev = o.external_state_variables.erase(pev);
      } else {
        ++pev;
      }
    }
    for (const auto& ev : i.external_state_variables) {
      update_field_holder(o.external_state_variables[ev.first], ev.second);
    }
    if (i.mass_density.has_value()) {
      if (!o.mass_density.has_value()) {
        o.mass_density = mgis::real{};
      }
      update_field_holder(*(o.mass_density), *(i.mass_density));
    } else {
      o.mass_density.reset();
    }
  }  // end of updateValues

  namespace internals {

    void extractScalarInternalStateVariable(
        mgis::span<mgis::real> o,
        const mgis::behaviour::MaterialStateManager& s,
        const mgis::size_type offset) {
      const auto stride = s.internal_state_variables_stride;
      auto* const p = o.data();
      const auto* const piv = s.internal_state_variables.data() + offset;
      for (mgis::size_type i = 0; i != s.n; ++i) {
        p[i] = piv[i * stride];
      }
    }  // end of extractScalarInternalStateVariable

    void extractInternalStateVariable(
        mgis::span<mgis::real> o,
        const mgis::behaviour::MaterialStateManager& s,
        const mgis::size_type nc,
        const mgis::size_type offset) {
      const auto stride = s.internal_state_variables_stride;
      auto* p = o.data();
      const auto* const piv = s.internal_state_variables.data() + offset;
      for (mgis::size_type i = 0; i != s.n; ++i) {
        const auto is = i * stride;
        for (mgis::size_type j = 0; j != nc; ++j, ++p) {
          *p = piv[is + j];
        }
      }
    }  // end of extractScalarInternalStateVariable

  }  // end of namespace internals

  void extractInternalStateVariable(
      mgis::span<mgis::real> o,
      const mgis::behaviour::MaterialStateManager& s,
      const std::string_view n) {
    const auto& iv = mgis::behaviour::getVariable(s.b.isvs, n);
    const auto nc = mgis::behaviour::getVariableSize(iv, s.b.hypothesis);
    const auto offset =
        mgis::behaviour::getVariableOffset(s.b.isvs, n, s.b.hypothesis);
    // checking compatibility
    if (o.size() != s.n * nc) {
      mgis::raise(
          "extractInternalStateVariable: "
          "unmatched number of integration points");
    }
    if (nc == 1) {
      mgis::behaviour::internals::extractScalarInternalStateVariable(o, s,
                                                                     offset);
    } else {
      mgis::behaviour::internals::extractInternalStateVariable(o, s, nc,
                                                               offset);
    }
  }  // end of extractInternalStateVariable

}  // end of namespace mgis::behaviour
