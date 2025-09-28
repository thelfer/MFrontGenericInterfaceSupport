/*!
 * \file   MGIS/Behaviour/MaterialFunctionManager.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/07/2025
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX
#define LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX

#include <iostream>

#include "MGIS/Behaviour/Variable.hxx"

namespace mgis::behaviour {

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  bool MaterialFunctionManager<SpaceType>::checkPreconditions(
      AbstractErrorHandler& eh,
      const std::shared_ptr<const SpaceType>& s,
      const Behaviour&) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid quadrature space");
    }
    return true;
  }  // end of checkPreconditions

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  MaterialFunctionManager<SpaceType>::MaterialFunctionManager(
      const std::shared_ptr<const SpaceType>& s, const Behaviour& behaviour)
      : MaterialFunctionManager(preconditions_check, s, behaviour) {}

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  template <bool doPreconditionsCheck>
  MaterialFunctionManager<SpaceType>::MaterialFunctionManager(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      const std::shared_ptr<const SpaceType>& s,
      const Behaviour& behaviour)
      : PreconditionsChecker<MaterialFunctionManager>(pcheck, s, behaviour),
        MaterialDataManager(behaviour,
                            static_cast<size_type>(getSpaceSize(*s))),
        qspace(s) {}  // end of MaterialFunctionManager

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  bool MaterialFunctionManager<SpaceType>::checkPreconditions(
      AbstractErrorHandler& eh,
      const std::shared_ptr<const SpaceType>& s,
      const std::shared_ptr<const Behaviour>& behaviour) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid quadrature space");
    }
    if (behaviour.get() == nullptr) {
      return eh.registerErrorMessage("invalid behaviour");
    }
    return true;
  }  // end of checkPreconditions

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  MaterialFunctionManager<SpaceType>::MaterialFunctionManager(
      const std::shared_ptr<const SpaceType>& s,
      const std::shared_ptr<const Behaviour>& behaviour)
      : MaterialFunctionManager(preconditions_check, s, behaviour) {}

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  template <bool doPreconditionsCheck>
  MaterialFunctionManager<SpaceType>::MaterialFunctionManager(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      const std::shared_ptr<const SpaceType>& s,
      const std::shared_ptr<const Behaviour>& behaviour)
      : PreconditionsChecker<MaterialFunctionManager>(pcheck, s, behaviour),
        MaterialDataManager(*(behaviour),
                            static_cast<size_type>(getSpaceSize(*s))),
        qspace(s),
        behaviour_ptr(behaviour) {}  // end of MaterialFunctionManager

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  const SpaceType& MaterialFunctionManager<SpaceType>::getSpace()
      const noexcept {
    return *(this->qspace);
  }  // end of getSpace

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::shared_ptr<const SpaceType>
  MaterialFunctionManager<SpaceType>::getSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getSpacePointer

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  mgis::function::SharedSpace<SpaceType>
  MaterialFunctionManager<SpaceType>::getSharedSpace() const noexcept {
    return {this->qspace};
  }  // end of getSpacePointer

  namespace internals {

    inline InvalidResult registerErrorMessage(AbstractErrorHandler& eh,
                                              std::string_view n,
                                              const char* const emsg1,
                                              const char* const emsg2) {
      auto* const ctx = dynamic_cast<Context*>(&eh);
      if (ctx != nullptr) {
        return ctx->registerErrorMessage(std::string(emsg1) + " '" +
                                         std::string{n} + "'");
      }
      return eh.registerErrorMessage(emsg2);
    }  // end of registerErrorMessage

    template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
    std::optional<mgis::function::FunctionView<
        mgis::function::SharedSpace<SpaceType>,
        fixed_size_dynamic_stride_data_layout_description<N>>>
    makeFunctionView(AbstractErrorHandler& eh,
                     const mgis::function::SharedSpace<SpaceType>& qspace,
                     const std::vector<Variable>& variables,
                     std::span<real> values,
                     std::string_view n,
                     const Hypothesis h,
                     const size_type stride,
                     const char* const emsg1,
                     const char* const emsg2) {
      using ::mgis::function::FunctionDataLayout;
      using ::mgis::function::FunctionView;
      if (!contains(variables, n)) {
        return registerErrorMessage(eh, n, emsg1, emsg2);
      }
      const auto& v = getVariable(variables, n);
      const auto vs = getVariableSize(v, h);
      const auto vo = getVariableOffset(variables, n, h);
      if constexpr (N != dynamic_extent) {
        if (N != vs) {
          return registerErrorMessage(eh, n,
                                      "invalid variable size for variable",
                                      "invalid variable size");
        }
      }
      if (getSpaceSize(qspace) == 0) {
        if constexpr (N == dynamic_extent) {
          return FunctionView<mgis::function::SharedSpace<SpaceType>>(
              qspace, {}, vs, stride);
        } else {
          return FunctionView<
              mgis::function::SharedSpace<SpaceType>,
              fixed_size_dynamic_stride_data_layout_description<N>>(qspace, {},
                                                                    stride);
        }
      }
      if constexpr (N == dynamic_extent) {
        return FunctionView<mgis::function::SharedSpace<SpaceType>>(
            qspace, values.subspan(vo), vs, stride);
      } else {
        return FunctionView<
            mgis::function::SharedSpace<SpaceType>,
            fixed_size_dynamic_stride_data_layout_description<N>>(
            qspace, values.subspan(vo), stride);
      }
    }  // end of makeFunctionView

    template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
    std::optional<mgis::function::FunctionView<
        mgis::function::SharedSpace<SpaceType>,
        fixed_size_dynamic_stride_data_layout_description<N>,
        false>>
    makeImmutableFunctionView(
        AbstractErrorHandler& eh,
        const mgis::function::SharedSpace<SpaceType>& qspace,
        const std::vector<Variable>& variables,
        std::span<const real> values,
        std::string_view n,
        const Hypothesis h,
        const size_type stride,
        const char* const emsg1,
        const char* const emsg2) {
      using ::mgis::function::FunctionDataLayout;
      using ::mgis::function::FunctionView;
      if (!contains(variables, n)) {
        return registerErrorMessage(eh, n, emsg1, emsg2);
      }
      const auto& v = getVariable(variables, n);
      const auto vs = getVariableSize(v, h);
      const auto vo = getVariableOffset(variables, n, h);
      if constexpr (N != dynamic_extent) {
        if (N != vs) {
          return registerErrorMessage(eh, n,
                                      "invalid variable size for variable",
                                      "invalid variable size");
        }
      }
      if (getSpaceSize(qspace) == 0) {
        if constexpr (N == dynamic_extent) {
          return FunctionView<mgis::function::SharedSpace<SpaceType>,
                              mgis::function::FunctionDataLayoutDescription{},
                              false>(qspace, {}, vs, stride);
        } else {
          return FunctionView<
              mgis::function::SharedSpace<SpaceType>,
              fixed_size_dynamic_stride_data_layout_description<N>, false>(
              qspace, {}, stride);
        }
      }
      if constexpr (N == dynamic_extent) {
        return FunctionView<mgis::function::SharedSpace<SpaceType>,
                            mgis::function::FunctionDataLayoutDescription{},
                            false>(qspace, values.subspan(vo), vs, stride);
      } else {
        return FunctionView<
            mgis::function::SharedSpace<SpaceType>,
            fixed_size_dynamic_stride_data_layout_description<N>, false>(
            qspace, values.subspan(vo), stride);
      }
    }  // end of makeFunctionView

  }  // end of namespace internals

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getGradient(AbstractErrorHandler& eh,
              MaterialFunctionManager<SpaceType>& m,
              std::string_view n,
              const TimeStepStage ts) {
    return getGradient<dynamic_extent>(eh, m, n, ts);
  }  // end of getGradient

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      mgis::function::FunctionDataLayoutDescription{},
      false>>
  getGradient(AbstractErrorHandler& eh,
              const MaterialFunctionManager<SpaceType>& m,
              std::string_view n,
              const TimeStepStage ts) {
    return getGradient<dynamic_extent>(eh, m, n, ts);
  }  // end of getGradient

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getGradient(AbstractErrorHandler& eh,
              MaterialFunctionManager<SpaceType>& m,
              std::string_view n,
              const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeFunctionView<N>(
        eh, m.getSharedSpace(), m.b.gradients, s.gradients, n, m.b.hypothesis,
        s.gradients_stride, "no gradient named",
        "no gradient with the given name");
  }  // end of getGradient

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getGradient(AbstractErrorHandler& eh,
              const MaterialFunctionManager<SpaceType>& m,
              std::string_view n,
              const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeImmutableFunctionView<N>(
        eh, m.getSharedSpace(), m.b.gradients, s.gradients, n, m.b.hypothesis,
        s.gradients_stride, "no gradient named",
        "no gradient with the given name");
  }  // end of getGradient

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getThermodynamicForce(AbstractErrorHandler& eh,
                        MaterialFunctionManager<SpaceType>& m,
                        std::string_view n,
                        const TimeStepStage ts) {
    return getThermodynamicForce<dynamic_extent>(eh, m, n, ts);
  }  // end of getThermodynamicForce

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      mgis::function::FunctionDataLayoutDescription{},
      false>>
  getThermodynamicForce(AbstractErrorHandler& eh,
                        const MaterialFunctionManager<SpaceType>& m,
                        std::string_view n,
                        const TimeStepStage ts) {
    return getThermodynamicForce<dynamic_extent>(eh, m, n, ts);
  }  // end of getThermodynamicForce

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getThermodynamicForce(AbstractErrorHandler& eh,
                        MaterialFunctionManager<SpaceType>& m,
                        std::string_view n,
                        const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeFunctionView<N>(
        eh, m.getSharedSpace(), m.b.thermodynamic_forces,
        s.thermodynamic_forces, n, m.b.hypothesis,
        s.thermodynamic_forces_stride, "no thermodynamic force named",
        "no thermodynamic force with the given name");
  }  // end of getThermodynamicForce

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getThermodynamicForce(AbstractErrorHandler& eh,
                        const MaterialFunctionManager<SpaceType>& m,
                        std::string_view n,
                        const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeImmutableFunctionView<N>(
        eh, m.getSharedSpace(), m.b.thermodynamic_forces,
        s.thermodynamic_forces, n, m.b.hypothesis,
        s.thermodynamic_forces_stride, "no thermodynamic force named",
        "no thermodynamic force with the given name");
  }  // end of getThermodynamicForce

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getInternalStateVariable(AbstractErrorHandler& eh,
                           MaterialFunctionManager<SpaceType>& m,
                           std::string_view n,
                           const TimeStepStage ts) {
    return getInternalStateVariable<dynamic_extent>(eh, m, n, ts);
  }  // end of getInternalStateVariable

  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      mgis::function::FunctionDataLayoutDescription{},
      false>>
  getInternalStateVariable(AbstractErrorHandler& eh,
                           const MaterialFunctionManager<SpaceType>& m,
                           std::string_view n,
                           const TimeStepStage ts) {
    return getInternalStateVariable<dynamic_extent>(eh, m, n, ts);
  }  // end of getInternalStateVariable

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getInternalStateVariable(AbstractErrorHandler& eh,
                           MaterialFunctionManager<SpaceType>& m,
                           std::string_view n,
                           const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeFunctionView<N>(
        eh, m.getSharedSpace(), m.b.isvs, s.internal_state_variables, n,
        m.b.hypothesis, s.internal_state_variables_stride,
        "no internal state variable named",
        "no internal state variable with the given name");
  }  // end of getInternalStateVariable

  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getInternalStateVariable(AbstractErrorHandler& eh,
                           const MaterialFunctionManager<SpaceType>& m,
                           std::string_view n,
                           const TimeStepStage ts) {
    auto& s = (ts == ets) ? m.s1 : m.s0;
    return internals::makeImmutableFunctionView<N>(
        eh, m.getSharedSpace(), m.b.isvs, s.internal_state_variables, n,
        m.b.hypothesis, s.internal_state_variables_stride,
        "no internal state variable named",
        "no internal state variable with the given name");
  }  // end of getInternalStateVariable

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX */
