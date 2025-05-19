/*!
 * \file   MGIS/Function/TensorView.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TENSORVIEW_IXX

#include <utility>

namespace mgis::function {

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  bool TensorView<FunctionType, TensorType>::checkPreconditions(
      const FunctionType& values) noexcept {
    return values.getNumberOfComponents() == compile_time_size<TensorType>;
  }  // end of checkPreconditions

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  TensorView<FunctionType, TensorType>::TensorView(FunctionType& values)
      : function(values) {}  // end of TensorView

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  bool TensorView<FunctionType, TensorType>::check(Context&) const noexcept {
    return checkPreconditions(this->function);
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  void TensorView<FunctionType, TensorType>::allocateWorkspace() {}

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  const typename TensorView<FunctionType, TensorType>::Space&
  TensorView<FunctionType, TensorType>::getSpace() const {
    return this->function.getSpace();
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr size_type
  TensorView<FunctionType, TensorType>::getNumberOfComponents() const noexcept {
    return compile_time_size<TensorType>;
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const element_index<Space>& i) const
      requires(internals::FunctionResultQuery<FunctionType>::b1) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(this->function.data(unsafe, i));
    } else {
      static_assert(std::same_as<decltype(static_cast<const FunctionType&>(
                                              this->function)(i)
                                              .data()),
                                 const real*>);
      return tfel::math::map<const TensorType>(this->function(i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires(internals::FunctionResultQuery<FunctionType>::b2) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, wk, i));
    } else {
      static_assert(
          std::same_as<decltype(this->function(wk, i).data()), const real*>);
      return tfel::math::map<const TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires(internals::FunctionResultQuery<FunctionType>::b3) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const cell_index<Space> re,
                 const quadrature_point_index<Space> ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, e, i));
    } else {
      static_assert(
          std::same_as<decltype(this->function(e, i).data()), const real*>);
      return tfel::math::map<const TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires(internals::FunctionResultQuery<FunctionType>::b4) {
    constexpr auto has_data_method = requires(
        const FunctionType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space> re, const quadrature_point_index<Space> ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, wk, e, i));
    } else {
      static_assert(
          std::same_as<decltype(this->function(wk, e, i).data()), const real*>);
      return tfel::math::map<const TensorType>(this->function(wk, e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const element_index<Space>& i)  //
      requires(internals::FunctionResultQuery<FunctionType>::b1) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, i));
    } else {
      static_assert(std::same_as<decltype(this->function(i).data()), real*>);
      return tfel::math::map<TensorType>(this->function(i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const element_workspace<Space>& wk,
      const element_index<Space>& i)  //
      requires(internals::FunctionResultQuery<FunctionType>::b2) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, i));
    } else {
      static_assert(
          std::same_as<decltype(this->function(wk, i).data()), real*>);
      return tfel::math::map<TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires(internals::FunctionResultQuery<FunctionType>::b3) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, e, i));
    } else {
      static_assert(std::same_as<decltype(this->function(e, i).data()), real*>);
      return tfel::math::map<TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto TensorView<FunctionType, TensorType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires(internals::FunctionResultQuery<FunctionType>::b4) {
    constexpr auto has_data_method = requires(
        FunctionType & rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, e, i));
    } else {
      static_assert(
          std::same_as<decltype(this->function(wk, e, i).data()), real*>);
      return tfel::math::map<TensorType>(this->function(wk, e, i).data());
    }
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORVIEW_IXX */
