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
#ifdef MGIS_HAVE_HDF5
#include "MGIS/Utilities/HDF5Support.hxx"
#endif /* MGIS_HAVE_HDF5 */
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis::behaviour {

  MaterialStateManager::MutableFieldHolder&
  MaterialStateManager::MutableFieldHolder::operator=(
      const mgis::real v) noexcept {
    this->value = v;
    return *this;
  }  // end of value

  MaterialStateManager::FieldHolder&
  MaterialStateManager::FieldHolder::operator=(const mgis::real v) noexcept {
    this->operator=(MutableFieldHolder{v});
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
                       std::span<mgis::real> evalues, const size_type vs,
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
      std::map<std::string, MaterialStateManager::FieldHolder, std::less<>>& m,
      std::string_view n) noexcept {
    return m[std::string{n}];
  }  // end of getFieldHolder

  static std::map<std::string, MaterialStateManager::FieldHolder, std::less<>>::
      const_iterator
      getFieldHolderIterator(const std::map<std::string,
                                            MaterialStateManager::FieldHolder,
                                            std::less<>>& m,
                             std::string_view n) {
    // #if __cplusplus > 201103L
    //       return m.find(n);
    // #else  /* __cplusplus > 201103L */
    return m.find(n);
    // #endif /* __cplusplus > 201103L */
  }  // end of getFieldHolder

  void setMaterialProperty(MaterialStateManager& m,
                           std::string_view n,
                           const real v,
                           const MaterialStateManager::UpdatePolicy p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    setMaterialProperty(ctx, m, n, v, p) | or_raise;
  }  // end of setMaterialProperty

  bool setMaterialProperty(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      const real v,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    const auto& mp = *omp;
    if (mp.type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "setMaterialProperty: "
          "invalid material property "
          "(only scalar material property is supported)");
    }
    auto& fh = getFieldHolder(m.material_properties, n);
    fh = MaterialStateManager::MutableFieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    return true;
  }  // end of setMaterialProperty

  void setMaterialProperty(MaterialStateManager& m,
                           std::string_view n,
                           std::span<real> v,
                           const MaterialStateManager::StorageMode s,
                           const MaterialStateManager::UpdatePolicy p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    setMaterialProperty(ctx, m, n, v, s, p) | or_raise;
  }  // end of setMaterialProperty

  bool setMaterialProperty(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      std::span<real> v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    const auto& mp = *omp;
    if (mp.type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "invalid material property (only scalar material property is "
          "supported)");
    }
    if (static_cast<mgis::size_type>(v.size()) != m.n) {
      return ctx.registerErrorMessage(
          "invalid number of values (does not match the number of integration "
          "points)");
    }
    auto& fh = getFieldHolder(m.material_properties, n);
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      return setMaterialProperty(ctx, m, n, std::span<const real>{v}, s, p);
    }
    fh = MaterialStateManager::MutableFieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    return true;
  }  // end of setMaterialProperty

  bool setMaterialProperty(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      std::span<const real> v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    const auto& mp = *omp;
    if (mp.type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "invalid material property "
          "(only scalar material property is supported)");
    }
    if (static_cast<mgis::size_type>(v.size()) != m.n) {
      return ctx.registerErrorMessage(
          "invalid number of values (does not match the number of integration "
          "points)");
    }
    auto& fh = getFieldHolder(m.material_properties, n);
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      fh = MaterialStateManager::MutableFieldHolder{
          .value = std::vector<real>{v.begin(), v.end()},
          .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    } else {
      if (p != MaterialStateManager::UPDATE) {
        return ctx.registerErrorMessage(
            "update policy must be equal to `NOUPDATE` when using a "
            "non-modifiable external storage");
      }
      fh = v;
    }
    return true;
  }  // end of setMaterialProperty

  bool unsetMaterialProperty(Context& ctx,
                             MaterialStateManager& m,
                             std::string_view n) noexcept {
    const auto omp = getVariable(ctx, m.b.mps, n);
    if (isInvalid(omp)) {
      return false;
    }
    auto p = m.material_properties.find(n);
    if (p != m.material_properties.end()) {
      m.material_properties.erase(p);
    }
    return true;
  }  // end of unsetMaterialProperty

  bool isMaterialPropertyDefined(const MaterialStateManager& m,
                                 std::string_view n) {
    const auto p = getFieldHolderIterator(m.material_properties, n);
    return p != m.material_properties.end();
  }  // end of isMaterialPropertyDefined

  bool isMaterialPropertyUniform(const MaterialStateManager& m,
                                 std::string_view n) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return isMaterialPropertyUniform(ctx, m, n) | or_raise;
  }  // end of isMaterialPropertyUniform

  std::optional<bool> isMaterialPropertyUniform(Context& ctx,
                                                const MaterialStateManager& m,
                                                std::string_view n) noexcept {
    using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
    const auto p = getFieldHolderIterator(m.material_properties, n);
    if (p == m.material_properties.end()) {
      return ctx.registerErrorMessage("no material property named '" +
                                      std::string{n} + "' defined");
    }
    if (!std::holds_alternative<MutableFieldHolder>(p->second)) {
      return false;
    }
    return std::holds_alternative<real>(
        std::get<MutableFieldHolder>(p->second).value);
  }  // end of isMaterialPropertyUniform

  void setMassDensity(MaterialStateManager& m,
                      const real v,
                      const MaterialStateManager::UpdatePolicy p) noexcept {
    m.mass_density = MaterialStateManager::MutableFieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
  }  // end of setMassDensity

  MGIS_EXPORT void setMassDensity(MaterialStateManager& m,
                                  std::span<real> v,
                                  const MaterialStateManager::StorageMode s,
                                  const MaterialStateManager::UpdatePolicy p) {
    mgis::raise_if(static_cast<mgis::size_type>(v.size()) != m.n,
                   "setMassDensity: invalid number of values "
                   "(does not match the number of integration points)");
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      m.mass_density = MaterialStateManager::MutableFieldHolder{
          .value = std::vector<real>{v.begin(), v.end()},
          .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    } else {
      m.mass_density = MaterialStateManager::MutableFieldHolder{
          .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    }
  }  // end of setMassDensity

  bool isMassDensityDefined(const MaterialStateManager& m) noexcept {
    return m.mass_density.has_value();
  }  // end of isMassDensityDefined

  bool isMassDensityUniform(const MaterialStateManager& m) {
    auto ctx = Context{};
    return isMassDensityUniform(ctx, m) | ctx.getThrowingFailureHandler();
  }  // end of isMassDensityUniform

  std::optional<bool> isMassDensityUniform(
      Context& ctx, const MaterialStateManager& m) noexcept {
    if (!isMassDensityDefined(m)) {
      return ctx.registerErrorMessage("the mass density is undefined");
    }
    return std::holds_alternative<real>(m.mass_density->value);
  }  // end of isMassDensityUniform

  void setExternalStateVariable(MaterialStateManager& m,
                                std::string_view n,
                                const real v,
                                const MaterialStateManager::UpdatePolicy p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    setExternalStateVariable(ctx, m, n, v, p) | or_raise;
  }  // end of setExternalStateVariable

  void setExternalStateVariable(MaterialStateManager& m,
                                std::string_view n,
                                std::span<real> v,
                                const MaterialStateManager::StorageMode s,
                                const MaterialStateManager::UpdatePolicy p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    setExternalStateVariable(ctx, m, n, v, s, p) | or_raise;
  }  // end of setExternalStateVariable

  bool setExternalStateVariable(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      const real v,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto oesv = getVariable(ctx, m.b.esvs, n);
    if (isInvalid(oesv)) {
      return false;
    }
    if (oesv->type != Variable::SCALAR) {
      return ctx.registerErrorMessage(
          "setExternalStateVariable: "
          "invalid external state variable "
          "(only scalar external state variable is supported)");
    }
    auto& fh = getFieldHolder(m.external_state_variables, n);
    fh = MaterialStateManager::MutableFieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    return true;
  }  // end of setExternalStateVariable

  bool setExternalStateVariable(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      std::span<real> v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto oesv = getVariable(ctx, m.b.esvs, n);
    if (isInvalid(oesv)) {
      return false;
    }
    const auto ovs = getVariableSize(ctx, *oesv, m.b.hypothesis);
    if (isInvalid(ovs)) {
      return false;
    }
    if (((static_cast<mgis::size_type>(v.size()) != m.n * (*ovs)) &&
         (static_cast<mgis::size_type>(v.size()) != (*ovs)))) {
      return ctx.registerErrorMessage(
          "setExternalStateVariable: invalid number of values");
    }
    auto& fh = getFieldHolder(m.external_state_variables, n);
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      return setExternalStateVariable(ctx, m, n, std::span<const real>{v}, s,
                                      p);
    }
    fh = MaterialStateManager::MutableFieldHolder{
        .value = v, .shall_be_updated = (p == MaterialStateManager::UPDATE)};
    return true;
  }  // end of setExternalStateVariable

  bool setExternalStateVariable(
      Context& ctx,
      MaterialStateManager& m,
      std::string_view n,
      std::span<const real> v,
      const MaterialStateManager::StorageMode s,
      const MaterialStateManager::UpdatePolicy p) noexcept {
    const auto oesv = getVariable(ctx, m.b.esvs, n);
    if (isInvalid(oesv)) {
      return false;
    }
    const auto ovs = getVariableSize(ctx, *oesv, m.b.hypothesis);
    if (isInvalid(ovs)) {
      return false;
    }
    if (((static_cast<mgis::size_type>(v.size()) != m.n * (*ovs)) &&
         (static_cast<mgis::size_type>(v.size()) != (*ovs)))) {
      return ctx.registerErrorMessage(
          "setExternalStateVariable: invalid number of values");
    }
    auto& fh = getFieldHolder(m.external_state_variables, n);
    if (s == MaterialStateManager::LOCAL_STORAGE) {
      if (v.size() == 1u) {
        fh = MaterialStateManager::MutableFieldHolder{
            .value = v[0],
            .shall_be_updated = (p == MaterialStateManager::UPDATE)};
      } else {
        fh = MaterialStateManager::MutableFieldHolder{
            .value = std::vector<real>{v.begin(), v.end()},
            .shall_be_updated = (p == MaterialStateManager::UPDATE)};
      }
    } else {
      if (p != MaterialStateManager::UPDATE) {
        return ctx.registerErrorMessage(
            "update policy must be equal to `NOUPDATE` when using a "
            "non-modifiable external storage");
      }
      fh = v;
    }
    return true;
  }  // end of setExternalStateVariable

  bool unsetExternalStateVariable(Context& ctx,
                                  MaterialStateManager& m,
                                  std::string_view n) noexcept {
    const auto oesv = getVariable(ctx, m.b.esvs, n);
    if (isInvalid(oesv)) {
      return false;
    }
    auto p = m.external_state_variables.find(n);
    if (p != m.external_state_variables.end()) {
      m.external_state_variables.erase(p);
    }
    return true;
  }  // end of unsetExternalStateVariable

  bool isExternalStateVariableDefined(const MaterialStateManager& m,
                                      std::string_view n) {
    const auto p = getFieldHolderIterator(m.external_state_variables, n);
    return p != m.external_state_variables.end();
  }  // end of isExternalStateVariableDefined

  bool isExternalStateVariableUniform(const MaterialStateManager& m,
                                      std::string_view n) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return isExternalStateVariableUniform(ctx, m, n) | or_raise;
  }  // end of isExternalStateVariableUniform

  std::optional<bool> isExternalStateVariableUniform(
      Context& ctx,
      const MaterialStateManager& m,
      std::string_view n) noexcept {
    using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
    const auto p = getFieldHolderIterator(m.external_state_variables, n);
    if (p == m.external_state_variables.end()) {
      return ctx.registerErrorMessage("no external state variable named '" +
                                      std::string{n} + "' defined");
    }
    if (!std::holds_alternative<MutableFieldHolder>(p->second)) {
      return false;
    }
    return std::holds_alternative<real>(
        std::get<MutableFieldHolder>(p->second).value);
  }  // end of isExternalStateVariableUniform

  [[nodiscard]] static bool updateMutableFieldHolderFromSpan(
      Context& ctx,
      MaterialStateManager& o,
      MaterialStateManager::MutableFieldHolder& dest,
      std::span<const real> values) noexcept {
    if (!dest.shall_be_updated) {
      return true;
    }
    auto& to = dest.value;
    if (std::holds_alternative<std::span<mgis::real>>(to)) {
      // reuse existing memory
      auto& to_v = std::get<std::span<mgis::real>>(to);
      if (values.size() != to_v.size()) {
        return ctx.registerErrorMessage("arrays' size does not match");
      }
      std::copy(values.begin(), values.end(), to_v.begin());
    } else if (std::holds_alternative<std::vector<mgis::real>>(to)) {
      // reuse existing memory
      auto& to_v = std::get<std::vector<mgis::real>>(to);
      if (to_v.size() != values.size()) {
        if (to_v.size() * o.n == values.size()) {
          to_v.resize(values.size());
        }
      }
      if (values.size() != to_v.size()) {
        return ctx.registerErrorMessage("arrays' size does not match");
      }
      std::copy(values.begin(), values.end(), to_v.begin());
    } else {
      // to contains a real value, so overwrite it with a new vector
      to = std::vector<mgis::real>(values.begin(), values.end());
    }
    return true;
  }

  [[nodiscard]] static bool updateMutableFieldHolder(
      Context& ctx,
      MaterialStateManager& o,
      MaterialStateManager::MutableFieldHolder& dest,
      const MaterialStateManager::MutableFieldHolder& src) noexcept {
    if (!dest.shall_be_updated) {
      return true;
    }
    const auto& from = src.value;
    if (std::holds_alternative<mgis::real>(from)) {
      dest.value = std::get<mgis::real>(from);
    } else if (std::holds_alternative<std::vector<mgis::real>>(from)) {
      return updateMutableFieldHolderFromSpan(
          ctx, o, dest, std::get<std::vector<mgis::real>>(from));
    } else {
      return updateMutableFieldHolderFromSpan(
          ctx, o, dest, std::get<std::span<mgis::real>>(from));
    }
    return true;
  }  // end of updateMutableFieldHolder

  [[nodiscard]] static bool updateFieldHolder(
      Context& ctx,
      MaterialStateManager& o,
      MaterialStateManager::FieldHolder& dest,
      const MaterialStateManager::FieldHolder& src) noexcept {
    using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
    if (std::holds_alternative<std::monostate>(dest)) {
      return ctx.registerErrorMessage("uninitialized field holder given");
    }
    if (std::holds_alternative<std::monostate>(src)) {
      return ctx.registerErrorMessage("uninitialized field holder given");
    }
    if (std::holds_alternative<std::span<const real>>(dest)) {
      return true;
    }
    auto& fh = std::get<MutableFieldHolder>(dest);
    if (!fh.shall_be_updated) {
      return true;
    }
    if (std::holds_alternative<MutableFieldHolder>(src)) {
      return updateMutableFieldHolder(ctx, o, fh,
                                      std::get<MutableFieldHolder>(src));
    }
    return updateMutableFieldHolderFromSpan(
        ctx, o, fh, std::get<std::span<const real>>(src));
  }  // end of updateFieldHolder

  [[nodiscard]] static bool updateFieldHolder(
      Context& ctx,
      MaterialStateManager& o,
      MaterialStateManager::FieldHolder& dest,
      const std::string& n,
      const MaterialStateManager::FieldHolder& src) noexcept {
    if (!updateFieldHolder(ctx, o, dest, src)) {
      return ctx.registerErrorMessage("updating field holder '" + n +
                                      "' failed");
    }
    return true;
  }  // end of updateFieldHolder

  bool updateValues(Context& ctx,
                    MaterialStateManager& o,
                    const MaterialStateManager& i) noexcept {
    using FieldHolder = MaterialStateManager::FieldHolder;
    using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
    auto update_span = [&ctx](
                           std::span<real>& to,
                           const std::span<const real>& from) noexcept -> bool {
      if (from.size() != to.size()) {
        return ctx.registerErrorMessage("arrays' size does not match");
      }
      std::copy(from.begin(), from.end(), to.begin());
      return true;
    };  // end update_span

    auto check_mps =
        [&ctx](
            const Behaviour& b,
            const std::map<std::string, FieldHolder, std::less<>>& mps) noexcept
        -> bool {
      for (const auto& mp : mps) {
        if (!contains(b.mps, mp.first)) {
          return ctx.registerErrorMessage(
              "material property '" + mp.first +
              "' defined in the material state manager is not defined "
              " by the behaviour");
        }
      }
      return true;
    };
    if (&i.b != &o.b) {
      return ctx.registerErrorMessage(
          "the material state managers do not holds the same behaviour");
    }
    if (!check_mps(o.b, i.material_properties)) {
      return false;
    }
    if (!check_mps(o.b, o.material_properties)) {
      return false;
    }
    if (!update_span(o.gradients, i.gradients)) {
      return false;
    }
    if (!update_span(o.thermodynamic_forces, i.thermodynamic_forces)) {
      return false;
    }
    if (!update_span(o.internal_state_variables, i.internal_state_variables)) {
      return false;
    }
    if (!update_span(o.stored_energies, i.stored_energies)) {
      return false;
    }
    if (!update_span(o.dissipated_energies, i.dissipated_energies)) {
      return false;
    }
    auto pmp = o.material_properties.begin();
    while (pmp != o.material_properties.end()) {
      if (i.material_properties.count(pmp->first) == 0) {
        pmp = o.material_properties.erase(pmp);
      } else {
        ++pmp;
      }
    }
    for (const auto& mp : i.material_properties) {
      if (o.material_properties.find(mp.first) == o.material_properties.end()) {
        if (std::holds_alternative<std::span<const real>>(mp.second)) {
          continue;
        }
        o.material_properties.insert({mp.first, MutableFieldHolder{real{}}});
      }
      if (!updateFieldHolder(ctx, o, o.material_properties.at(mp.first),
                             mp.first, mp.second)) {
        return false;
      }
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
      if (o.external_state_variables.find(ev.first) ==
          o.external_state_variables.end()) {
        if (std::holds_alternative<std::span<const real>>(ev.second)) {
          continue;
        }
        o.external_state_variables.insert(
            {ev.first, MutableFieldHolder{real{}}});
      }
      if (!updateFieldHolder(ctx, o, o.external_state_variables.at(ev.first),
                             ev.first, ev.second)) {
        return false;
      }
    }
    if (i.mass_density.has_value()) {
      if (!o.mass_density.has_value()) {
        o.mass_density = MutableFieldHolder{.value = mgis::real{}};
      }
      if (!updateMutableFieldHolder(ctx, o, *(o.mass_density),
                                    *(i.mass_density))) {
        return false;
      }
    } else {
      o.mass_density.reset();
    }
    return true;
    }  // end of updateValues

    void updateValues(MaterialStateManager & o, const MaterialStateManager& i) {
      auto ctx = Context{};
      updateValues(ctx, o, i) | ctx.getThrowingFailureHandler();
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
        std::span<mgis::real> o, const mgis::behaviour::MaterialStateManager& s,
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
      const MaterialStateManager::MutableFieldHolder& f,
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

  [[nodiscard]] static bool save(
      Context& ctx,
      H5::Group& g,
      const std::string& n,
      const MaterialStateManager::FieldHolder& f,
      const MaterialStateManagerSavingOptions& opts) noexcept {
    using namespace mgis::utilities::hdf5;
    if (std::holds_alternative<std::monostate>(f)) {
      return ctx.registerErrorMessage("uninitialized field holder '" + n + "'");
    }
    if (std::holds_alternative<std::span<const real>>(f)) {
      return write(ctx, g, n, std::get<std::span<const real>>(f),
                   opts.allow_overwrite);
    }
    return save(ctx, g, n,
                std::get<MaterialStateManager::MutableFieldHolder>(f), opts);
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

  std::optional<MaterialStateManagerRestoreOptions>
  getGreedyMaterialStateManagerRestoreOptions(Context& ctx,
                                              const Behaviour& b,
                                              const H5::Group& g) noexcept {
    using namespace mgis::utilities::hdf5;
    const auto odatasets = getDataSetNames(ctx, g);
    if (isInvalid(odatasets)) {
      return {};
    }
    auto contains = [odatasets](std::string_view n) noexcept {
      return std::find(odatasets->begin(), odatasets->end(), n) !=
             odatasets->end();
    };
    auto select_ignored_variables = [&g, &ctx](
                                        const std::vector<Variable>& variables,
                                        const std::string& gn) {
      auto og = openGroup(ctx, g, gn);
      auto ignored_variables = std::vector<std::string>{};
      if (isInvalid(og)) {
        return ignored_variables;
      }
      const auto oldatasets = getDataSetNames(ctx, *og);
      if (isInvalid(oldatasets)) {
        return ignored_variables;
      }
      for (const auto& v : variables) {
        const auto found = std::find(oldatasets->begin(), oldatasets->end(),
                                     v.name) != oldatasets->end();
        if (!found) {
          ignored_variables.push_back(v.name);
        }
      }
      return ignored_variables;
    };
    return MaterialStateManagerRestoreOptions{
        .restore_gradients = contains("gradients"),
        .restore_thermodynamic_forces = contains("thermodynamic_forces"),
        .restore_stored_energies = contains("stored_energies"),
        .restore_dissipated_energies = contains("dissipated_energies"),
        .restore_internal_state_variables =
            contains("internal_state_variables"),
        .restore_mass_densities = contains("mass_density"),
        .restore_material_properties = subGroupExists(g, "material_properties"),
        .ignored_material_properties =
            select_ignored_variables(b.mps, "material_properties"),
        .restore_external_state_variables =
            subGroupExists(g, "external_state_variables"),
        .ignored_external_state_variables =
            select_ignored_variables(b.esvs, "external_state_variables")};
  }  // end of getGreedyMaterialStateManagerRestoreOptions

  [[nodiscard]] static bool restoreScalarMutableFieldHolder(
      Context& ctx,
      MaterialStateManager::MutableFieldHolder& f,
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
  }  // end of restoreScalarMutableFieldHolder

  [[nodiscard]] static bool restoreMutableFieldHolder(
      Context& ctx,
      MaterialStateManager::MutableFieldHolder& f,
      const H5::Group& g,
      const std::string& n,
      const size_type ng,
      const size_type vsize) noexcept {
    using namespace mgis::utilities::hdf5;
    if (vsize == 1) {
      return restoreScalarMutableFieldHolder(ctx, f, g, n, ng);
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
  }  // end of restoreMutableFieldHolder

  bool restore(Context& ctx,
               MaterialStateManager& s,
               const H5::Group& g,
               const MaterialStateManagerRestoreOptions& opts) noexcept {
    using namespace mgis::utilities::hdf5;
    using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
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
      if (!restoreScalarMutableFieldHolder(ctx, *(s.mass_density), g,
                                           "mass_density", s.n)) {
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
        if (std::find(opts.ignored_material_properties.begin(),
                      opts.ignored_material_properties.end(),
                      mp.name) != opts.ignored_material_properties.end()) {
          continue;
        }
        const auto ovsize = getVariableSize(ctx, mp, s.b.hypothesis);
        if (isInvalid(ovsize)) {
          return ctx.registerErrorMessage("restoring material property '" +
                                          mp.name + "' failed");
        }
        auto& f = s.material_properties[mp.name];
        if (std::holds_alternative<std::monostate>(f)) {
          return ctx.registerErrorMessage("uninitialized material property '" +
                                          mp.name + "'");
        }
        if (std::holds_alternative<MutableFieldHolder>(f)) {
          auto& mfh = std::get<MutableFieldHolder>(f);
          if (!restoreMutableFieldHolder(ctx, mfh, *ogmp, mp.name, s.n,
                                         *ovsize)) {
            return ctx.registerErrorMessage("restoring material property '" +
                                            mp.name + "' failed");
          }
        } else {
          if (!opts.ignore_constant_variables) {
            return ctx.registerErrorMessage("material property '" + mp.name +
                                            "' can't be restored: it has been "
                                            "declared as a constant value");
          }
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
        if (std::find(opts.ignored_external_state_variables.begin(),
                      opts.ignored_external_state_variables.end(), esv.name) !=
            opts.ignored_external_state_variables.end()) {
          continue;
        }
        const auto ovsize = getVariableSize(ctx, esv, s.b.hypothesis);
        if (isInvalid(ovsize)) {
          return ctx.registerErrorMessage(
              "restoring external state variable '" + esv.name + "' failed");
        }
        auto& f = s.external_state_variables[esv.name];
        if (std::holds_alternative<std::monostate>(f)) {
          return ctx.registerErrorMessage(
              "uninitialized external state variable'" + esv.name + "'");
        }
        if (std::holds_alternative<MutableFieldHolder>(f)) {
          auto& mfh = std::get<MutableFieldHolder>(f);
          if (!restoreMutableFieldHolder(ctx, mfh, *ogesv, esv.name, s.n,
                                         *ovsize)) {
            return ctx.registerErrorMessage(
                "restoring external state variable '" + esv.name + "' failed");
          }
        } else {
          if (!opts.ignore_constant_variables) {
            return ctx.registerErrorMessage("external state variable '" +
                                            esv.name +
                                            "' can't be restored: it has been "
                                            "declared as a constant value");
          }
        }
      }
    }
    return true;
  }  // end of restore

#endif /* MGIS_HAVE_HDF5 */

  }  // end of namespace mgis::behaviour
