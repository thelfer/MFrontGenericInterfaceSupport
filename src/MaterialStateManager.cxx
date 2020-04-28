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

namespace mgis {

  namespace behaviour {

    MaterialStateManager::MaterialStateManager(const Behaviour& behaviour,
                                               const size_type s)
        : gradients_stride(
              getArraySize(behaviour.gradients, behaviour.hypothesis)),
          thermodynamic_forces_stride(getArraySize(
              behaviour.thermodynamic_forces, behaviour.hypothesis)),
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
      init(this->thermodynamic_forces, this->thermodynamic_forces_values,
           this->thermodynamic_forces_stride);
      init(this->internal_state_variables,
           this->internal_state_variables_values,
           this->internal_state_variables_stride);
      init(this->stored_energies, this->stored_energies_values, 1u);
      init(this->dissipated_energies, this->dissipated_energies_values, 1u);
    }  // end of MaterialStateManager::MaterialStateManager

    MaterialStateManager::MaterialStateManager(
        const Behaviour& behaviour,
        const size_type s,
        const MaterialStateManagerInitializer& i)
        : gradients_stride(
              getArraySize(behaviour.gradients, behaviour.hypothesis)),
          thermodynamic_forces_stride(getArraySize(
              behaviour.thermodynamic_forces, behaviour.hypothesis)),
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
      init(this->thermodynamic_forces, this->thermodynamic_forces_values,
           i.thermodynamic_forces, this->thermodynamic_forces_stride,
           "thermodynamic forces");
      init(this->internal_state_variables,
           this->internal_state_variables_values, i.internal_state_variables,
           this->internal_state_variables_stride, "internal state variables");
      init(this->stored_energies, this->stored_energies_values,
           i.stored_energies, 1u, "stored energies");
      init(this->dissipated_energies, this->dissipated_energies_values,
           i.dissipated_energies, 1u, "dissipated energies");
    }  // end of MaterialStateManager::MaterialStateManager

    MaterialStateManager::~MaterialStateManager() = default;

    static MaterialStateManager::FieldHolder& getFieldHolder(
        std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
// #if __cplusplus > 201103L
//       return m[n];
// #else  /* __cplusplus > 201103L */
      return m[n.to_string()];
// #endif /* __cplusplus > 201103L */
    } // end of getFieldHolder

    static std::map<std::string, MaterialStateManager::FieldHolder>::iterator
    getFieldHolderIterator(
        std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
// #if __cplusplus > 201103L
//       return m.find(n);
// #else  /* __cplusplus > 201103L */
      return m.find(n.to_string());
      // #endif /* __cplusplus > 201103L */
    } // end of getFieldHolder

    static std::map<std::string,
                    MaterialStateManager::FieldHolder>::const_iterator
    getFieldHolderIterator(
        const std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
// #if __cplusplus > 201103L
//       return m.find(n);
// #else  /* __cplusplus > 201103L */
      return m.find(n.to_string());
// #endif /* __cplusplus > 201103L */
    } // end of getFieldHolder

    void setMaterialProperty(MaterialStateManager& m,
                             const mgis::string_view& n,
                             const real v) {
      const auto mp = getVariable(m.b.mps, n);
      mgis::raise_if(mp.type != Variable::SCALAR,
                     "setMaterialProperty: "
                     "invalid material property "
                     "(only scalar material property is supported)");
      getFieldHolder(m.material_properties,n)= v;
    }  // end of setMaterialProperty

    MGIS_EXPORT void setMaterialProperty(
        MaterialStateManager& m,
        const mgis::string_view& n,
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
        getFieldHolder(m.material_properties,n)= std::vector<real>{v.begin(), v.end()};
      } else {
        getFieldHolder(m.material_properties,n)= v;
      }
    }  // end of setMaterialProperty

    bool isMaterialPropertyDefined(const MaterialStateManager& m,
                                        const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.material_properties,n);
      return p != m.material_properties.end();
    }  // end of isMaterialPropertyDefined

    bool isMaterialPropertyUniform(const MaterialStateManager& m,
                                   const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.material_properties, n);
      if (p == m.material_properties.end()) {
        mgis::raise(
            "isMaterialPropertyUniform: "
            "no material property named '" +
            n.to_string() + "' defined");
      }
      return mgis::holds_alternative<real>(p->second);
    }  // end of isMaterialPropertyUniform

    real& getUniformMaterialProperty(MaterialStateManager& m,
                                     const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.material_properties, n);
      mgis::raise_if(p == m.material_properties.end(),
                     "getUniformMaterialProperty: "
                     "no material property named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(!mgis::holds_alternative<real>(p->second),
                     "getUniformMaterialProperty: "
                     "material property '" +
                         n.to_string() + "' is not uniform");
      return mgis::get<real>(p->second);
    }  // end of getUniformMaterialProperty

    const real& getUniformMaterialProperty(const MaterialStateManager& m,
                                           const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.material_properties, n);
      mgis::raise_if(p == m.material_properties.end(),
                     "getUniformMaterialProperty: "
                     "no material property named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(!mgis::holds_alternative<real>(p->second),
                     "getUniformMaterialProperty: "
                     "material property '" +
                         n.to_string() + "' is not uniform");
      return mgis::get<real>(p->second);
    }  // end of getUniformMaterialProperty

    mgis::span<real> getNonUniformMaterialProperty(MaterialStateManager& m,
                                                   const mgis::string_view& n) {
      using index_type = span<real>::index_type;
      const auto p = getFieldHolderIterator(m.material_properties, n);
      mgis::raise_if(p == m.material_properties.end(),
                     "getNonUniformMaterialProperty: "
                     "no material property named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(mgis::holds_alternative<real>(p->second),
                     "getNonUniformMaterialProperty: "
                     "material property '" +
                         n.to_string() + "' is uniform");
      if (mgis::holds_alternative<std::vector<real>>(p->second)) {
        auto& values = mgis::get<std::vector<real>>(p->second);
        return {&values[0], static_cast<index_type>(values.size())};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getNonUniformMaterialProperty

    mgis::span<const real> getNonUniformMaterialProperty(
        const MaterialStateManager& m, const mgis::string_view& n) {
      using index_type = span<const real>::index_type;
      const auto p = getFieldHolderIterator(m.material_properties, n);
      mgis::raise_if(p == m.material_properties.end(),
                     "getNonUniformMaterialProperty: "
                     "no material property named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(mgis::holds_alternative<real>(p->second),
                     "getNonUniformMaterialProperty: "
                     "material property '" +
                         n.to_string() + "' is uniform");
      if (mgis::holds_alternative<std::vector<real>>(p->second)) {
        const auto& values = mgis::get<std::vector<real>>(p->second);
        return {&values[0], static_cast<index_type>(values.size())};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getNonUniformMaterialProperty

    void setExternalStateVariable(MaterialStateManager& m,
                                  const mgis::string_view& n,
                                  const real v) {
      const auto esv = getVariable(m.b.esvs, n);
      mgis::raise_if(esv.type != Variable::SCALAR,
                     "setExternalStateVariable: "
                     "invalid external state variable "
                     "(only scalar external state variable is supported)");
      getFieldHolder(m.external_state_variables,n) = v;
    }  // end of setExternalStateVariable

    MGIS_EXPORT void setExternalStateVariable(
        MaterialStateManager& m,
        const mgis::string_view& n,
        const mgis::span<real>& v,
        const MaterialStateManager::StorageMode s) {
      const auto esv = getVariable(m.b.esvs, n);
      mgis::raise_if(esv.type != Variable::SCALAR,
                     "setExternalStateVariable: "
                     "invalid external state variable "
                     "(only scalar external state variable is supported)");
      mgis::raise_if(static_cast<mgis::size_type>(v.size()) != m.n,
                     "setExternalStateVariable: invalid number of values "
                     "(does not match the number of integration points)");
      if (s == MaterialStateManager::LOCAL_STORAGE) {
        getFieldHolder(m.external_state_variables, n) =
            std::vector<real>{v.begin(), v.end()};
      } else {
        getFieldHolder(m.external_state_variables,n) = v;
      }
    }  // end of setExternalStateVariable

    bool isExternalStateVariableDefined(const MaterialStateManager& m,
                                        const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      return p != m.external_state_variables.end();
    }  // end of isExternalStateVariableDefined

    bool isExternalStateVariableUniform(const MaterialStateManager& m,
                                   const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      mgis::raise_if(p == m.external_state_variables.end(),
                     "isExternalStateVariableUniform: "
                     "no external state variable named '" +
                         n.to_string() + "' defined");
      return mgis::holds_alternative<real>(p->second);
    }  // end of isExternalStateVariableUniform

    real& getUniformExternalStateVariable(MaterialStateManager& m,
                                          const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      mgis::raise_if(p == m.external_state_variables.end(),
                     "getUniformExternalStateVariable: "
                     "no external state variable named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(!mgis::holds_alternative<real>(p->second),
                     "getUniformExternalStateVariable: "
                     "external state variable '" +
                         n.to_string() + "' is not uniform");
      return mgis::get<real>(p->second);
    }  // end of getUniformExternalStateVariable

    const real& getUniformExternalStateVariable(const MaterialStateManager& m,
                                                const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      mgis::raise_if(p == m.external_state_variables.end(),
                     "getUniformExternalStateVariable: "
                     "no external state variable named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(!mgis::holds_alternative<real>(p->second),
                     "getUniformExternalStateVariable: "
                     "external state variable '" +
                         n.to_string() + "' is not uniform");
      return mgis::get<real>(p->second);
    }  // end of getUniformExternalStateVariable

    mgis::span<real> getNonUniformExternalStateVariable(
        MaterialStateManager& m, const mgis::string_view& n) {
      using index_type = span<real>::index_type;
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      mgis::raise_if(p == m.external_state_variables.end(),
                     "getNonUniformExternalStateVariable: "
                     "no external state variable named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(mgis::holds_alternative<real>(p->second),
                     "getNonUniformExternalStateVariable: "
                     "external state variable '" +
                         n.to_string() + "' is uniform");
      if (mgis::holds_alternative<std::vector<real>>(p->second)) {
        auto& values = mgis::get<std::vector<real>>(p->second);
        return {&values[0],static_cast<index_type>(values.size())};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getUniformExternalStateVariable

    mgis::span<const real> getNonUniformExternalStateVariable(
        const MaterialStateManager& m, const mgis::string_view& n) {
      const auto p = getFieldHolderIterator(m.external_state_variables,n);
      using index_type = span<const real>::index_type;
      mgis::raise_if(p == m.external_state_variables.end(),
                     "getNonUniformExternalStateVariable: "
                     "no external state variable named '" +
                         n.to_string() + "' defined");
      mgis::raise_if(mgis::holds_alternative<real>(p->second),
                     "getNonUniformExternalStateVariable: "
                     "external state variable '" +
                         n.to_string() + "' is uniform");
      if (mgis::holds_alternative<std::vector<real>>(p->second)) {
        const auto& values = mgis::get<std::vector<real>>(p->second);
        return {&values[0], static_cast<index_type>(values.size())};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getUniformExternalStateVariable

    void update_values(MaterialStateManager& o, const MaterialStateManager& i) {
      auto check_size = [](const mgis::size_type s1, const mgis::size_type s2) {
        if (s1 != s2) {
          mgis::raise(
              "mgis::behaviour::update_values: "
              "arrays' size does not match");
        }
      };  // end of check_size
      auto update_span = [&check_size](mgis::span<real>& to,
                                       const mgis::span<const real>& from) {
        check_size(from.size(), to.size());
        std::copy(from.begin(), from.end(), to.begin());
      };  // end update_span
      auto update_field_holder = [&check_size](
          MaterialStateManager::FieldHolder& to,
          const MaterialStateManager::FieldHolder& from) {
        if (mgis::holds_alternative<mgis::real>(from)) {
          to = mgis::get<mgis::real>(from);
        } else if(mgis::holds_alternative<std::vector<mgis::real>>(from)){
          const auto& from_v = mgis::get<std::vector<mgis::real>>(from);
          if (mgis::holds_alternative<mgis::span<mgis::real>>(to)) {
            // reuse existing memory
            auto& to_v = mgis::get<mgis::span<mgis::real>>(to);
            check_size(from_v.size(), to_v.size());
            std::copy(from_v.begin(), from_v.end(), to_v.begin());
          } else if(mgis::holds_alternative<std::vector<mgis::real>>(to)) {
            // reuse existing memory
            auto& to_v = mgis::get<std::vector<mgis::real>>(to);
            check_size(from_v.size(), to_v.size());
            std::copy(from_v.begin(), from_v.end(), to_v.begin());
          } else {
            // to contains a real value, so overwrite it with a new vector
            to = mgis::get<std::vector<mgis::real>>(from);
          }
        } else {
          const auto from_v = mgis::get<mgis::span<mgis::real>>(from);
          if (mgis::holds_alternative<mgis::span<mgis::real>>(to)) {
            // reuse existing memory
            auto to_v = mgis::get<mgis::span<mgis::real>>(to);
            check_size(from_v.size(), to_v.size());
            std::copy(from_v.begin(), from_v.end(),
                      to_v.begin());
          } else if (mgis::holds_alternative<std::vector<mgis::real>>(to)) {
            // reuse existing memory
            auto to_v = mgis::get<std::vector<mgis::real>>(to);
            check_size(from_v.size(), to_v.size());
            std::copy(from_v.begin(), from_v.end(),
                      to_v.begin());
          } else {
            to = from_v;
          }
        }
      };  // end of update_field_holder
      auto check_mps = [](
          const Behaviour& b,
          const std::map<std::string, MaterialStateManager::FieldHolder>& mps) {
        for (const auto& mp : mps) {
          auto find_mp = [&mp](const Variable& d) {
            return mp.first == d.name;
          };
          if (std::find_if(b.mps.begin(), b.mps.end(), find_mp) ==
              b.mps.end()) {
            mgis::raise(
                "mgis::behaviour::update_values: "
                "material property '" +
                mp.first +
                "' defined in the material state manager is not defined "
                " by the behaviour");
          }
        }
      };
      if (&i.b != &o.b) {
        mgis::raise(
            "mgis::behaviour::update_values: the material state managers "
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
    }  // end of update_values

  }  // end of namespace behaviour

}  // end of namespace mgis
