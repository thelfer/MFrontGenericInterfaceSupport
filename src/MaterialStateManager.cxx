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

#include <utility>
#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis::behaviour {

  MaterialStateManager::FieldHolder&
  MaterialStateManager::FieldHolder::operator=(const mgis::real v) noexcept {
    this->value = v;
    return *this;
  }  // end of value

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
    auto init = [this](std::span<mgis::real>& view,
                       std::vector<mgis::real>& values, const size_type vs) {
      constexpr const auto zero = real{0};
      values.resize(this->n * vs, zero);
      view = std::span<mgis::real>(values);
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
    auto init = [this](std::span<mgis::real>& view,
                       std::vector<mgis::real>& values,
                       const std::span<mgis::real>& evalues, const size_type vs,
                       const char* const vn) {
      if (evalues.empty()) {
        constexpr const auto zero = real{0};
        values.resize(this->n * vs, zero);
        view = std::span<mgis::real>(values);
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

  [[nodiscard]] static MaterialStateManager::FieldHolder& getFieldHolder(
      std::map<std::string, MaterialStateManager::FieldHolder>& m,
      const std::string_view& n) noexcept {
    // #if __cplusplus > 201103L
    //       return m[n];
    // #else  /* __cplusplus > 201103L */
    return m[std::string{n}];
    // #endif /* __cplusplus > 201103L */
  }  // end of getFieldHolder

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
                           const real v,
                           const MaterialStateManager::UpdatePolicy p) {
    const auto mp = getVariable(m.b.mps, n);
    mgis::raise_if(mp.type != Variable::SCALAR,
                   "setMaterialProperty: "
                   "invalid material property "
                   "(only scalar material property is supported)");
    getFieldHolder(m.material_properties,
                   n) = MaterialStateManager::FieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
  }  // end of setMaterialProperty

  bool setMaterialProperty(
      Context& ctx,
      MaterialStateManager& m,
      const std::string_view& n,
      const real v,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    const auto& mp = *(*omp);
    if (mp.type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "setMaterialProperty: "
          "invalid material property "
          "(only scalar material property is supported)");
    }
    getFieldHolder(m.material_properties,
                   n) = MaterialStateManager::FieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    return true;
  }  // end of setMaterialProperty

  void setMaterialProperty(MaterialStateManager& m,
                           const std::string_view& n,
                           const std::span<real>& v,
                           const MaterialStateManager::StorageMode s,
                           const MaterialStateManager::UpdatePolicy p) {
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
          MaterialStateManager ::FieldHolder{
              .value = std::vector<real>{v.begin(), v.end()},
              .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    } else {
      getFieldHolder(m.material_properties,
                     n) = MaterialStateManager ::FieldHolder{
          .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    }
  }  // end of setMaterialProperty

  bool setMaterialProperty(
      Context& ctx,
      MaterialStateManager& m,
      const std::string_view& n,
      const std::span<real>& v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    const auto& mp = *(*omp);
    if (mp.type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "setMaterialProperty: "
          "invalid material property "
          "(only scalar material property is supported)");
    }
    if (static_cast<mgis::size_type>(v.size()) != m.n) {
      return ctx.registerErrorMessage(
          "setMaterialProperty: invalid number of values "
          "(does not match the number of integration points)");
    }
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      getFieldHolder(m.material_properties, n) =
          MaterialStateManager ::FieldHolder{
              .value = std::vector<real>{v.begin(), v.end()},
              .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    } else {
      getFieldHolder(m.material_properties,
                     n) = MaterialStateManager ::FieldHolder{
          .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    }
    return true;
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
    return std::holds_alternative<real>(p->second.value);
  }  // end of isMaterialPropertyUniform

  void setMassDensity(MaterialStateManager& m,
                      const real v,
                      const MaterialStateManager::UpdatePolicy p) {
    m.mass_density = MaterialStateManager::FieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
  }  // end of setMassDensity

  MGIS_EXPORT void setMassDensity(MaterialStateManager& m,
                                  const std::span<real>& v,
                                  const MaterialStateManager::StorageMode s,
                                  const MaterialStateManager::UpdatePolicy p) {
    mgis::raise_if(static_cast<mgis::size_type>(v.size()) != m.n,
                   "setMassDensity: invalid number of values "
                   "(does not match the number of integration points)");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      m.mass_density = MaterialStateManager::FieldHolder{
          .value = std::vector<real>{v.begin(), v.end()},
          .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    } else {
      m.mass_density = MaterialStateManager::FieldHolder{
          .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    }
  }  // end of setMassDensity

  bool isMassDensityDefined(const MaterialStateManager& m) {
    return m.mass_density.has_value();
  }  // end of isMassDensityDefined

  bool isMassDensityUniform(const MaterialStateManager& m) {
    if (!isMassDensityDefined(m)) {
      mgis::raise("isMassDensityUniform: the mass density is undefined");
    }
    return std::holds_alternative<real>(m.mass_density->value);
  }  // end of isMassDensityUniform

  void setExternalStateVariable(MaterialStateManager& m,
                                const std::string_view& n,
                                const real v,
                                const MaterialStateManager::UpdatePolicy p) {
    const auto esv = getVariable(m.b.esvs, n);
    mgis::raise_if(esv.type != Variable::SCALAR,
                   "setExternalStateVariable: "
                   "invalid external state variable "
                   "(only scalar external state variable is supported)");
    getFieldHolder(m.external_state_variables,
                   n) = MaterialStateManager::FieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
  }  // end of setExternalStateVariable

  MGIS_EXPORT void setExternalStateVariable(
      MaterialStateManager& m,
      const std::string_view& n,
      const std::span<real>& v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) {
    const auto esv = getVariable(m.b.esvs, n);
    const auto vs = getVariableSize(esv, m.b.hypothesis);
    mgis::raise_if(((static_cast<mgis::size_type>(v.size()) != m.n * vs) &&
                    (static_cast<mgis::size_type>(v.size()) != vs)),
                   "setExternalStateVariable: invalid number of values");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      if (v.size() == 1u) {
        getFieldHolder(m.external_state_variables, n) =
            MaterialStateManager::FieldHolder{
                .value = v[0],
                .shall_be_updated = (p == MaterialStateManager::UPDATE)};
      } else {
        getFieldHolder(m.external_state_variables, n) =
            MaterialStateManager::FieldHolder{
                .value = std::vector<real>{v.begin(), v.end()},
                .shall_be_updated = (p == MaterialStateManager::UPDATE)};
      }
    } else {
      getFieldHolder(m.external_state_variables,
                     n) = MaterialStateManager::FieldHolder{
          .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
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
    return std::holds_alternative<real>(p->second.value);
  }  // end of isExternalStateVariableUniform

  void updateValues(MaterialStateManager& o, const MaterialStateManager& i) {
    auto check_size = [](const mgis::size_type s1, const mgis::size_type s2) {
      if (s1 != s2) {
        mgis::raise(
            "mgis::behaviour::updateValues: "
            "arrays' size does not match");
      }
    };  // end of check_size
    auto update_span = [&check_size](std::span<real>& to,
                                     const std::span<const real>& from) {
      check_size(from.size(), to.size());
      std::copy(from.begin(), from.end(), to.begin());
    };  // end update_span
    auto update_field_holder =
        [&o, &check_size](MaterialStateManager::FieldHolder& dest,
                          const MaterialStateManager::FieldHolder& src) {
          if (!dest.shall_be_updated) {
            return;
          }
          auto& to = dest.value;
          const auto& from = src.value;
          if (std::holds_alternative<mgis::real>(from)) {
            to = std::get<mgis::real>(from);
          } else if (std::holds_alternative<std::vector<mgis::real>>(from)) {
            const auto& from_v = std::get<std::vector<mgis::real>>(from);
            if (std::holds_alternative<std::span<mgis::real>>(to)) {
              // reuse existing memory
              auto& to_v = std::get<std::span<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else if (std::holds_alternative<std::vector<mgis::real>>(to)) {
              // reuse existing memory
              auto& to_v = std::get<std::vector<mgis::real>>(to);
              if (to_v.size() != from_v.size()) {
                if (to_v.size() * o.n == from_v.size()) {
                  to_v.resize(from_v.size());
                }
              }
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else {
              // to contains a real value, so overwrite it with a new vector
              to = std::get<std::vector<mgis::real>>(from);
            }
          } else {
            const auto from_v = std::get<std::span<mgis::real>>(from);
            if (std::holds_alternative<std::span<mgis::real>>(to)) {
              // reuse existing memory
              auto to_v = std::get<std::span<mgis::real>>(to);
              check_size(from_v.size(), to_v.size());
              std::copy(from_v.begin(), from_v.end(), to_v.begin());
            } else if (std::holds_alternative<std::vector<mgis::real>>(to)) {
              // reuse existing memory
              auto& to_v = std::get<std::vector<mgis::real>>(to);
              if (to_v.size() != from_v.size()) {
                if (to_v.size() * o.n == from_v.size()) {
                  to_v.resize(from_v.size());
                }
              }
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
        o.mass_density =
            MaterialStateManager::FieldHolder{.value = mgis::real{}};
      }
      update_field_holder(*(o.mass_density), *(i.mass_density));
    } else {
      o.mass_density.reset();
    }
  }  // end of updateValues

  namespace internals {
    static void extractScalarInternalStateVariable(
        std::span<mgis::real> o,
        const mgis::behaviour::MaterialStateManager& s,
        const mgis::size_type offset) {
      const auto stride = s.internal_state_variables_stride;
      auto* const p = o.data();
      const auto* const piv = s.internal_state_variables.data() + offset;
      for (mgis::size_type i = 0; i != s.n; ++i) {
        p[i] = piv[i * stride];
      }
    }  // end of extractScalarInternalStateVariable

    static void extractInternalStateVariable(
        std::span<mgis::real> o,
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
      std::span<mgis::real> o,
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

#ifdef MGIS_HAVE_HDF5

  [[nodiscard]] static bool save(
      Context& ctx,
      H5::Group& g,
      const std::string& n,
      const MaterialStateManager::FieldHolder& f,
      const MaterialStateManagerSavingOptions& opts) noexcept {
    using namespace mgis::utilities::hdf5;
    if (std::holds_alternative<real>(f.value)) {
      return write(ctx, g, n, std::get<real>(f.value), opts.allow_overwrite);
    } else if (std::holds_alternative<std::span<real>>(f.value)) {
      return write(ctx, g, n, std::get<std::span<real>>(f.value),
                   opts.allow_overwrite);
    }
    return write(ctx, g, n, std::get<std::vector<real>>(f.value),
                 opts.allow_overwrite);
  }  // end of save

  bool save(Context& ctx,
            H5::Group& g,
            const MaterialStateManager& s,
            const MaterialStateManagerSavingOptions& opts) noexcept {
    using namespace mgis::utilities::hdf5;
    if (opts.save_gradients) {
      if (!write(ctx, g, "gradients", s.gradients, opts.allow_overwrite)) {
        return false;
      }
    }
    if (opts.save_thermodynamic_forces) {
      if (!write(ctx, g, "thermodynamic_forces", s.thermodynamic_forces,
                 opts.allow_overwrite)) {
        return false;
      }
    }
    if (s.mass_density.has_value()) {
      if (!save(ctx, g, "mass_density", *(s.mass_density), opts)) {
        return false;
      }
    }
    //
    if (opts.save_material_properties) {
      auto og_mp = createGroup(ctx, g, "material_properties");
      if (isInvalid(og_mp)) {
        return false;
      }
      for (const auto& [n, mp] : s.material_properties) {
        if (!save(ctx, *og_mp, n, mp, opts)) {
          return false;
        }
      }
    }
    //
    if (!s.b.isvs.empty()) {
      if (!write(ctx, g, "internal_state_variables", s.internal_state_variables,
                 opts.allow_overwrite)) {
        return false;
      }
    }
    if ((opts.save_stored_energies) && (!s.stored_energies.empty())) {
      if (!write(ctx, g, "stored_energies", s.stored_energies,
                 opts.allow_overwrite)) {
        return false;
      }
    }
    if ((opts.save_dissipated_energies) && (!s.dissipated_energies.empty())) {
      if (!write(ctx, g, "dissipated_energies", s.dissipated_energies,
                 opts.allow_overwrite)) {
        return false;
      }
    }
    //
    if (opts.save_external_state_variables) {
      auto og_esv = createGroup(ctx, g, "external_state_variables");
      if (isInvalid(og_esv)) {
        return false;
      }
      for (const auto& [n, esv] : s.external_state_variables) {
        if (!save(ctx, *og_esv, n, esv, opts)) {
          return false;
        }
      }
    }
    return true;
  }  // end of save

  [[nodiscard]] static bool restoreScalarFieldHolder(
      Context& ctx,
      MaterialStateManager::FieldHolder& f,
      const H5::Group& g,
      const std::string& n,
      const size_type ng) noexcept {
    using namespace mgis::utilities::hdf5;
    auto odataset = openDataSet(ctx, g, n);
    if (isInvalid(odataset)) {
      return false;
    }
    try {
      const auto s = odataset->getSpace();
      if (s.getSimpleExtentNdims() != 1) {
        return ctx.registerErrorMessage("unexpected multidimensional array");
      }
      hsize_t dims[1];
      s.getSimpleExtentDims(dims);
      if (dims[0] == 1) {
        // uniform value
        auto value = real{};
        odataset->read(&value, getNativeType<real>());
        f = value;
        return true;
      }
      if (dims[0] != ng) {
        return ctx.registerErrorMessage(
            "unexpected data size for '" + n + "' (expected '" +
            std::to_string(ng) + "', got '" + std::to_string(dims[0]) + "')");
      }
      if (std::holds_alternative<std::span<real>>(f.value)) {
        auto& values = std::get<std::span<real>>(f.value);
        const auto fsize = values.size();
        if (fsize != ng) {
          return ctx.registerErrorMessage(
              "invalid storage for '" + n + "' (expected '" +
              std::to_string(ng) + "', got '" + std::to_string(fsize) + "')");
        }
        odataset->read(values.data(), getNativeType<real>());
      } else if (std::holds_alternative<std::vector<real>>(f.value)) {
        auto& values = std::get<std::vector<real>>(f.value);
        if (values.empty()) {
          values.resize(ng);
        }
        const auto fsize = values.size();
        if (fsize != ng) {
          return ctx.registerErrorMessage(
              "invalid storage for '" + n + "' (expected '" +
              std::to_string(ng) + "', got '" + std::to_string(fsize) + "')");
        }
        odataset->read(values.data(), getNativeType<real>());
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of restoreScalarFieldHolder

  [[nodiscard]] static bool restoreFieldHolder(
      Context& ctx,
      MaterialStateManager::FieldHolder& f,
      const H5::Group& g,
      const std::string& n,
      const size_type ng,
      const size_type vsize) noexcept {
    using namespace mgis::utilities::hdf5;
    if (vsize == 1) {
      return restoreScalarFieldHolder(ctx, f, g, n, ng);
    }
    if (std::holds_alternative<real>(f.value)) {
      // by default, f can be initialized to a real value
      f.value = std::vector<real>{};
    }
    auto odataset = openDataSet(ctx, g, n);
    if (isInvalid(odataset)) {
      return false;
    }
    try {
      const auto s = odataset->getSpace();
      if (s.getSimpleExtentNdims() != 1) {
        return ctx.registerErrorMessage("unexpected multidimensional array");
      }
      hsize_t dims[1];
      s.getSimpleExtentDims(dims);
      if (dims[0] == vsize) {
        // uniform value
        if (std::holds_alternative<std::span<real>>(f.value)) {
          auto& values = std::get<std::span<real>>(f.value);
          if (values.size() == vsize) {
            odataset->read(values.data(), getNativeType<real>());
            return true;
          }
          if (values.size() != ng * vsize) {
            return ctx.registerErrorMessage(
                "unexpected data size for '" + n + "' (expected '" +
                std::to_string(ng * vsize) + "', got '" +
                std::to_string(dims[0]) + "')");
          }
          // duplicate the uniform values
          auto tmp = std::vector<real>(vsize);
          odataset->read(tmp.data(), getNativeType<real>());
          for (size_type idx = 0; idx != ng; ++idx) {
            std::copy(tmp.begin(), tmp.end(), values.begin() + idx * vsize);
          }
          return true;
        }
        auto& values = std::get<std::vector<real>>(f.value);
        if (values.empty()) {
          values.resize(vsize);
        }
        odataset->read(values.data(), getNativeType<real>());
        return true;
      }
      if (dims[0] != vsize * ng) {
        return ctx.registerErrorMessage(
            "unexpected data size for '" + n + "' (expected '" +
            std::to_string(ng * vsize) + "', got '" + std::to_string(dims[0]) +
            "')");
      }
      if (std::holds_alternative<std::span<real>>(f.value)) {
        auto& values = std::get<std::span<real>>(f.value);
        const auto fsize = values.size();
        if (fsize != ng) {
          return ctx.registerErrorMessage(
              "invalid storage for '" + n + "' (expected '" +
              std::to_string(ng) + "', got '" + std::to_string(fsize) + "')");
        }
        odataset->read(values.data(), getNativeType<real>());
      } else if (std::holds_alternative<std::vector<real>>(f.value)) {
        auto& values = std::get<std::vector<real>>(f.value);
        if (values.empty()) {
          values.resize(ng);
        }
        const auto fsize = values.size();
        if (fsize != ng) {
          return ctx.registerErrorMessage(
              "invalid storage for '" + n + "' (expected '" +
              std::to_string(ng) + "', got '" + std::to_string(fsize) + "')");
        }
        odataset->read(values.data(), getNativeType<real>());
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of restoreFieldHolder

  bool restore(Context& ctx,
               MaterialStateManager& s,
               const H5::Group& g,
               const MaterialStateManagerRestoreOptions& opts) noexcept {
    using namespace mgis::utilities::hdf5;
    if (opts.restore_gradients) {
      if (!read(ctx, s.gradients, g, "gradients")) {
        return ctx.registerErrorMessage("restoring gradients failed");
      }
    }
    if (opts.restore_thermodynamic_forces) {
      if (!read(ctx, s.thermodynamic_forces, g, "thermodynamic_forces")) {
        return ctx.registerErrorMessage(
            "restoring thermodynamic forces failed");
      }
    }
    if (opts.restore_mass_densities) {
      if (!s.mass_density.has_value()) {
        s.mass_density = real{};
      }
      if (!restoreScalarFieldHolder(ctx, *(s.mass_density), g, "mass_density",
                                    s.n)) {
        return ctx.registerErrorMessage("restoring mass density failed");
      }
    }
    if (opts.restore_internal_state_variables) {
      if (!read(ctx, s.internal_state_variables, g,
                "internal_state_variables")) {
        return ctx.registerErrorMessage(
            "restoring internal state variables failed");
      }
    }
    if ((opts.restore_stored_energies) && (s.b.computesStoredEnergy)) {
      if (!read(ctx, s.stored_energies, g, "stored_energies")) {
        return ctx.registerErrorMessage("restoring stored energies failed");
        return false;
      }
    }
    if ((opts.restore_dissipated_energies) && (s.b.computesDissipatedEnergy)) {
      if (!read(ctx, s.dissipated_energies, g, "dissipated_energies")) {
        return ctx.registerErrorMessage("restoring dissipated energies failed");
      }
    }
    if ((opts.restore_material_properties) && (!s.b.mps.empty())) {
      if (!subGroupExists(g, "material_properties")) {
        return ctx.registerErrorMessage("no material properties saved");
      }
      auto ogmp = openGroup(ctx, g, "material_properties");
      if (isInvalid(ogmp)) {
        return false;
      }
      for (const auto& mp : s.b.mps) {
        const auto ovsize = getVariableSize(ctx, mp, s.b.hypothesis);
        if (isInvalid(ovsize)) {
          return ctx.registerErrorMessage("restoring material property '" +
                                          mp.name + "' failed");
        }
        auto& f = s.material_properties[mp.name];
        if (!restoreFieldHolder(ctx, f, *ogmp, mp.name, s.n, *ovsize)) {
          return ctx.registerErrorMessage("restoring material property '" +
                                          mp.name + "' failed");
        }
      }
    }
    if ((opts.restore_external_state_variables) && (!s.b.esvs.empty())) {
      if (!subGroupExists(g, "external_state_variables")) {
        return ctx.registerErrorMessage("no external state variables saved");
      }
      auto ogesv = openGroup(ctx, g, "external_state_variables");
      if (isInvalid(ogesv)) {
        return false;
      }
      for (const auto& esv : s.b.esvs) {
        const auto ovsize = getVariableSize(ctx, esv, s.b.hypothesis);
        if (isInvalid(ovsize)) {
          return ctx.registerErrorMessage(
              "restoring external state variable '" + esv.name + "' failed");
        }
        auto& f = s.external_state_variables[esv.name];
        if (!restoreFieldHolder(ctx, f, *ogesv, esv.name, s.n, *ovsize)) {
          return ctx.registerErrorMessage(
              "restoring external state variable '" + esv.name + "' failed");
        }
      }
    }
    return true;
  }  // end of restore

#endif /* MGIS_HAVE_HDF5 */

}  // end of namespace mgis::behaviour
