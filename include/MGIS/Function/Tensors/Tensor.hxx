/*!
 * \file   MGIS/Function/Tensors/Tensor.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use tensor evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TENSOR_HXX
#define LIB_MGIS_FUNCTION_TENSOR_HXX

#include <type_traits>
#include "TFEL/Math/fsarray.hxx"
#include "TFEL/Math/tvector.hxx"
#include "TFEL/Math/tmatrix.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"

namespace mgis::function::internals {

  //! \brief partial specialization for symmetric tensors
  template <unsigned short N>
  requires((N == 1) && (N == 2) && (N == 3))  //
      struct CompileTimeSize<tfel::math::stensor<N, real>> {
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>;
  };
  //! \brief partial specialization for tensors
  template <unsigned short N>
  requires((N == 1) && (N == 2) && (N == 3))  //
      struct CompileTimeSize<tfel::math::tensor<N, real>> {
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>;
  };

  //! \brief partial specialization for finite size array
  template <unsigned short N>
  struct CompileTimeSize<tfel::math::fsarray<N, real>> {
    static constexpr size_type value = N;
  };

  //! \brief partial specialization for tiny vectors
  template <unsigned short N>
  struct CompileTimeSize<tfel::math::tvector<N, real>> {
    static constexpr size_type value = N;
  };

  //! \brief partial specialization for tiny matrices
  template <unsigned short N, unsigned short M>
  struct CompileTimeSize<tfel::math::tmatrix<N, M, real>> {
    static constexpr size_type value = N * M;
  };

  //! \brief partial specialization for mutable and immutable views
  template <typename T>
  struct CompileTimeSize<tfel::math::View<T>>
      : CompileTimeSize<std::remove_cv_t<T>> {};

} // end of mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSOR_HXX */
