/*!
 * \file   FixedSizeModifier.ixx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_IXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) bool FixedSizeModifier<EvaluatorType, N>::checkPreconditions(
      const EvaluatorType& values) noexcept {
    return values.getNumberOfComponents() == N;
  }  // end of checkPreconditions

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) FixedSizeModifier<EvaluatorType, N>::FixedSizeModifier(
      const EvaluatorType& e)
      : evaluator(e) {}  // end of FixedSizeModifier

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) bool FixedSizeModifier<EvaluatorType, N>::check(
      Context&) const noexcept {
    return checkPreconditions(this->evaluator);
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0)  //
      void FixedSizeModifier<EvaluatorType, N>::allocateWorkspace() {
    return this->evaluator.allocateWorkspace();
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0)                                                 //
      const typename FixedSizeModifier<EvaluatorType, N>::Space&  //
      FixedSizeModifier<EvaluatorType, N>::getSpace() const {
    return internals::disambiguateGetSpace(this->evaluator);
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) constexpr size_type
      FixedSizeModifier<EvaluatorType, N>::getNumberOfComponents()
          const noexcept {
    return N;
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) auto FixedSizeModifier<EvaluatorType, N>::operator()(
      const element_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b1) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->evaluator.data(unsafe, i));
      } else {
        return *(this->evaluator(i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->evaluator.data(unsafe, i), N);
      } else {
        return std::span<const real, N>(this->evaluator(i).data(), N);
      }
    }
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) auto FixedSizeModifier<EvaluatorType, N>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b2) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->evaluator.data(unsafe, wk, i));
      } else {
        return *(this->evaluator(wk, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->evaluator.data(unsafe, wk, i), N);
      } else {
        return std::span<const real, N>(this->evaluator(wk, i).data(), N);
      }
    }
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) auto FixedSizeModifier<EvaluatorType, N>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b3) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->evaluator.data(unsafe, e, i));
      } else {
        return *(this->evaluator(e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->evaluator.data(unsafe, e, i), N);
      } else {
        return std::span<const real, N>(this->evaluator(e, i).data(), N);
      }
    }
  }

  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) auto FixedSizeModifier<EvaluatorType, N>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires(internals::EvaluatorResultQuery<EvaluatorType>::b4) {
    constexpr auto has_data_method = requires(
        const EvaluatorType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->evaluator.data(unsafe, wk, e, i));
      } else {
        return *(this->evaluator(wk, e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->evaluator.data(unsafe, wk, e, i),
                                        N);
      } else {
        return std::span<const real, N>(this->evaluator(wk, e, i).data(), N);
      }
    }
  }

  template <size_type N, typename EvaluatorType>
  auto view(const EvaluatorType& f) requires(
      (N > 0) && (N != dynamic_extent) &&  //
      (EvaluatorConcept<std::decay_t<EvaluatorType>>)) {
    return FixedSizeModifier<std::decay_t<EvaluatorType>, N>(f);
  }  // end of view

  template <EvaluatorConcept EvaluatorType, size_type N>
  const auto& getSpace(const FixedSizeModifier<EvaluatorType, N>& e){
    return e.getSpace();
  } // end of getSpace

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER!_IXX */
