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

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis {

  namespace behaviour {

    MaterialStateManager::MaterialStateManager(MaterialStateManager&&) =
        default;
    MaterialStateManager::MaterialStateManager(const MaterialStateManager&) =
        default;

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
      auto init = [this](std::vector<real>& values, const size_type vs) {
        constexpr const auto zero = real{0};
        values.resize(this->n * vs, zero);
      };
      init(this->gradients, this->gradients_stride);
      init(this->thermodynamic_forces, this->thermodynamic_forces_stride);
      init(this->internal_state_variables,
           this->internal_state_variables_stride);
      init(this->stored_energies, 1u);
      init(this->dissipated_energies, 1u);
    }  // end of MaterialStateManager::MaterialStateManager

    MaterialStateManager& MaterialStateManager::operator=(
        MaterialStateManager&& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialStateManager::operator=: unmatched behaviour");
      mgis::raise_if(src.n != this->n,
                     "MaterialStateManager::operator=: unmatched number "
                     "of integration points");
      if (&src != this) {
        this->gradients = std::move(src.gradients);
        this->thermodynamic_forces = std::move(src.thermodynamic_forces);
        this->material_properties = std::move(src.material_properties);
        this->internal_state_variables =
            std::move(src.internal_state_variables);
        this->stored_energies = std::move(src.stored_energies);
        this->dissipated_energies = std::move(src.dissipated_energies);
        this->external_state_variables =
            std::move(src.external_state_variables);
      }
      return *this;
    }  // end of MaterialStateManager::operator=

    MaterialStateManager& MaterialStateManager::operator=(
        const MaterialStateManager& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialStateManager::operator=: unmatched behaviour");
      mgis::raise_if(src.n != this->n,
                     "MaterialStateManager::operator=: unmatched number "
                     "of integration points");
      if (&src != this) {
        this->gradients = src.gradients;
        this->thermodynamic_forces = src.thermodynamic_forces;
        this->material_properties = src.material_properties;
        this->internal_state_variables = src.internal_state_variables;
        this->stored_energies = src.stored_energies;
        this->dissipated_energies = src.dissipated_energies;
        this->external_state_variables = src.external_state_variables;
      }
      return *this;
    }  // end of MaterialStateManager::operator=

    MaterialStateManager::~MaterialStateManager() = default;

    static MaterialStateManager::FieldHolder& getFieldHolder(
        std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
#if __cplusplus > 201103L
      return m[n];
#else  /* __cplusplus > 201103L */
      return m[n.to_string()];
#endif /* __cplusplus > 201103L */
    } // end of getFieldHolder

    static std::map<std::string, MaterialStateManager::FieldHolder>::iterator
    getFieldHolderIterator(
        std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
#if __cplusplus > 201103L
      return m.find(n);
#else  /* __cplusplus > 201103L */
      return m.find(n.to_string());
#endif /* __cplusplus > 201103L */
    } // end of getFieldHolder

    static std::map<std::string,
                    MaterialStateManager::FieldHolder>::const_iterator
    getFieldHolderIterator(
        const std::map<std::string, MaterialStateManager::FieldHolder>& m,
        const mgis::string_view& n) {
#if __cplusplus > 201103L
      return m.find(n);
#else  /* __cplusplus > 201103L */
      return m.find(n.to_string());
#endif /* __cplusplus > 201103L */
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
      mgis::raise_if(v.size() != m.n,
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
        return {&values[0],values.size()};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getNonUniformMaterialProperty

    mgis::span<const real> getNonUniformMaterialProperty(
        const MaterialStateManager& m, const mgis::string_view& n) {
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
        return {&values[0],values.size()};
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
      mgis::raise_if(v.size() != m.n,
                     "setExternalStateVariable: invalid number of values "
                     "(does not match the number of integration points)");
      if (s == MaterialStateManager::LOCAL_STORAGE) {
        getFieldHolder(m.external_state_variables,n) = std::vector<real>{v.begin(), v.end()};
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
        return {&values[0],values.size()};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getUniformExternalStateVariable

    mgis::span<const real> getNonUniformExternalStateVariable(
        const MaterialStateManager& m, const mgis::string_view& n) {
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
        const auto& values = mgis::get<std::vector<real>>(p->second);
        return {&values[0],values.size()};
      }
      return mgis::get<span<real>>(p->second);
    }  // end of getUniformExternalStateVariable

  }  // end of namespace behaviour

}  // end of namespace mgis
