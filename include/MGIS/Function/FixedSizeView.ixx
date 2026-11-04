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
          ConstructorArgumentType values)
      : FixedSizeView(preconditions_check, values) {}  // end of FixedSizeView

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      template <bool doPreconditionsCheck>
      constexpr FixedSizeView<FunctionType, N>::FixedSizeView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          ConstructorArgumentType values)
      : PreconditionsChecker<FixedSizeView>(pcheck, values),
        function(make_view(values)) {}  // end of FixedSizeView

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr bool FixedSizeView<FunctionType, N>::check(
          AbstractErrorHandler& eh) const {
    if (internals::disambiguateGetNumberOfComponents(this->function) != N) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of check

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr decltype(auto)
          FixedSizeView<FunctionType, N>::getSpace() const {
    return internals::disambiguateGetSpace(this->function);
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr size_type
      FixedSizeView<FunctionType, N>::getNumberOfComponents() const MGIS_NOEXCEPT {
    return N;
  }

  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0)  //
      constexpr auto FixedSizeView<FunctionType, N>::operator()(
          const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        const function_view<FunctionType>& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, i));
      } else {
        return *(getDataPointer(this->function(i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, i), N);
      } else {
        return std::span<const real, N>(getDataPointer(this->function(i)), N);
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
    constexpr auto has_data_method = requires(
        const function_view<FunctionType>& rf,
        const element_workspace<Space>& rwk, const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, i));
      } else {
        return *(getDataPointer(this->function(wk, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, wk, i), N);
      } else {
        return std::span<const real, N>(getDataPointer(this->function(wk, i)),
                                        N);
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
    constexpr auto has_data_method = requires(
        const function_view<FunctionType>& rf, const cell_index<Space>& re,
        const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, e, i));
      } else {
        return *(getDataPointer(this->function(e, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, e, i), N);
      } else {
        return std::span<const real, N>(getDataPointer(this->function(e, i)),
                                        N);
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
        const function_view<FunctionType>& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, e, i));
      } else {
        return *(getDataPointer(this->function(wk, e, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<const real, N>(this->function.data(unsafe, wk, e, i),
                                        N);
      } else {
        return std::span<const real, N>(
            getDataPointer(this->function(wk, e, i)), N);
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
    constexpr auto has_data_method = requires(
        function_view<FunctionType>& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, i));
      } else {
        return *(getDataPointer(this->function(i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, i), N);
      } else {
        return std::span<real, N>(getDataPointer(this->function(i)), N);
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
    constexpr auto has_data_method = requires(
        function_view<FunctionType>& rf, const element_workspace<Space>& rwk,
        const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, i));
      } else {
        return *(getDataPointer(this->function(wk, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, wk, i), N);
      } else {
        return std::span<real, N>(getDataPointer(this->function(wk, i)), N);
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
    constexpr auto has_data_method = requires(
        function_view<FunctionType>& rf, const cell_index<Space>& re,
        const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, e, i));
      } else {
        return *(getDataPointer(this->function(e, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, e, i), N);
      } else {
        return std::span<real, N>(getDataPointer(this->function(e, i)), N);
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
        function_view<FunctionType>& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (N == 1) {
      if constexpr (has_data_method) {
        return *(this->function.data(unsafe, wk, e, i));
      } else {
        return *(getDataPointer(this->function(wk, e, i)));
      }
    } else {
      if constexpr (has_data_method) {
        return std::span<real, N>(this->function.data(unsafe, wk, e, i), N);
      } else {
        return std::span<real, N>(getDataPointer(this->function(wk, e, i)), N);
      }
    }
  }

  template <size_type N, typename FunctionType>
  constexpr auto view(FunctionType& f) requires(
      (N > 0) && (N != dynamic_extent) &&              //
      (FunctionConcept<std::decay_t<FunctionType>>)&&  //
      ((!std::is_const_v<FunctionType>) ||
       (LightweightViewConcept<std::decay_t<FunctionType>>))) {
    return FixedSizeView<std::decay_t<FunctionType>, N>(f);
  }  // end of view

  template <size_type N, typename FunctionType>
  constexpr auto view(FunctionType& f) requires(
      (N > 0) && (N != dynamic_extent) &&                //
      (FunctionConcept<std::decay_t<FunctionType>>)&&    //
      (!EvaluatorConcept<std::decay_t<FunctionType>>)&&  //
      (!LightweightViewConcept<std::decay_t<FunctionType>>)&&  //
      (std::is_const_v<FunctionType>)&&                  //
      (requires(FunctionType & rf) { view(rf); })) {
    return view<N>(view(f));
  }  // end of view

  template <FunctionConcept FunctionType, size_type N>
  constexpr decltype(auto) getSpace(const FixedSizeView<FunctionType, N>& v) {
    return v.getSpace();
  }  // end of getSpace

  template <FunctionConcept FunctionType, size_type N>
  constexpr bool check(AbstractErrorHandler& eh,
                       const FixedSizeView<FunctionType, N>& v) {
    return v.check(eh);
  }  // end of check

  template <FunctionConcept FunctionType, size_type N>
  constexpr mgis::size_type getNumberOfComponents(
      const FixedSizeView<FunctionType, N>& v) MGIS_NOEXCEPT {
    return v.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX */
