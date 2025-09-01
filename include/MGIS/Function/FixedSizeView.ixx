/*!
 * \file   FixedSizeView.ixx
 * \brief
 * \author th202608
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr bool FixedSizeView<FunctionType, N>::checkPreconditions(
          AbstractErrorHandler& eh, const FunctionType& values) {
    if (internals::disambiguateGetNumberOfComponents(values) != N) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr FixedSizeView<FunctionType, N>::FixedSizeView(
          FunctionType& values)
      : FixedSizeView(preconditions_check, values) {}  // end of FixedSizeView

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr FixedSizeView<FunctionType, N>::FixedSizeView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          FunctionType& values)
      : PreconditionsChecker<FixedSizeView>(pcheck, values),
        function(values) {}  // end of FixedSizeView

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr bool FixedSizeView<FunctionType, N>::check(
          AbstractErrorHandler& ctx) const {
    return checkPreconditions(ctx, this->function);
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr decltype(auto)
          FixedSizeView<FunctionType, N>::getSpace() const {
    return internals::disambiguateGetSpace(this->function);
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr size_type
      FixedSizeView<FunctionType, N>::getNumberOfComponents() const noexcept {
    return N;
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr auto FixedSizeView<FunctionType, N>::operator()(
          const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, i));
      } else {
        return *(this->function(i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, i), N);
      } else {
        return std::span<const real, N>(this->function(i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr auto FixedSizeView<FunctionType, N>::operator()(
          const element_workspace<Space>& wk,
          const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b2) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, i));
      } else {
        return *(this->function(wk, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, wk, i), N);
      } else {
        return std::span<const real, N>(this->function(wk, i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr auto FixedSizeView<FunctionType, N>::operator()(
          const cell_index<Space>& e,
          const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b3) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, e, i));
      } else {
        return *(this->function(e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, e, i), N);
      } else {
        return std::span<const real, N>(this->function(e, i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr auto FixedSizeView<FunctionType, N>::operator()(
          const cell_workspace<Space>& wk,
          const cell_index<Space>& e,
          const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        const FunctionType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, e, i));
      } else {
        return *(this->function(wk, e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, wk, e, i),
                                        N);
      } else {
        return std::span<const real, N>(this->function(wk, e, i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr typename FixedSizeView<FunctionType, N>::mutable_value_type
      FixedSizeView<FunctionType, N>::operator()(
          const element_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, i));
      } else {
        return *(this->function(i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, i), N);
      } else {
        return std::span<real, N>(this->function(i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr typename FixedSizeView<FunctionType, N>::mutable_value_type
      FixedSizeView<FunctionType, N>::operator()(
          const element_workspace<Space>& wk,
          const element_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b2) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, i));
      } else {
        return *(this->function(wk, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, wk, i), N);
      } else {
        return std::span<real, N>(this->function(wk, i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr typename FixedSizeView<FunctionType, N>::mutable_value_type
      FixedSizeView<FunctionType, N>::operator()(
          const cell_index<Space>& e,
          const quadrature_point_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b3) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, e, i));
      } else {
        return *(this->function(e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, e, i), N);
      } else {
        return std::span<real, N>(this->function(e, i).data(), N);
      }
    }
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr typename FixedSizeView<FunctionType, N>::mutable_value_type
      FixedSizeView<FunctionType, N>::operator()(
          const cell_workspace<Space>& wk,
          const cell_index<Space>& e,
          const quadrature_point_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        FunctionType & rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, e, i));
      } else {
        return *(this->function(wk, e, i).data());
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, wk, e, i), N);
      } else {
        return std::span<real, N>(this->function(wk, e, i).data(), N);
      }
    }
  }

  template <size_type N, typename FunctionType>
  constexpr auto view(FunctionType& f) requires(
      (N > 0) && (N != dynamic_extent) &&              //
      (FunctionConcept<std::decay_t<FunctionType>>)&&  //
      (!std::is_rvalue_reference_v<FunctionType>)) {
    return FixedSizeView<std::decay_t<FunctionType>, N>(f);
  }  // end of view

  template <FunctionConcept FunctionType, size_type N>
  constexpr decltype(auto) getSpace(const FixedSizeView<FunctionType, N>& v) {
    return v.getSpace();
  }  // end of getSpace

  template <FunctionConcept FunctionType, size_type N>
  constexpr void allocateWorkspace(FixedSizeView<FunctionType, N>&) noexcept {
  }  // end of allocateWorkspace

  template <FunctionConcept FunctionType, size_type N>
  constexpr mgis::size_type getNumberOfComponents(
      const FixedSizeView<FunctionType, N>& v) noexcept {
    return v.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX */
