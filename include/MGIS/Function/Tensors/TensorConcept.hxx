/*!
 * \file   MGIS/Function/Tensors/TensorConcept.hxx
 * \brief
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use tensor evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TENSORCONCEPT_HXX
#define LIB_MGIS_FUNCTION_TENSORCONCEPT_HXX

#include <type_traits>
#include "TFEL/Math/fsarray.hxx"
#include "TFEL/Math/tvector.hxx"
#include "TFEL/Math/tmatrix.hxx"
#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/t2tot2.hxx"
#include "TFEL/Math/t2tost2.hxx"
#include "TFEL/Math/st2tot2.hxx"
#include "TFEL/Math/st2tost2.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"

namespace mgis::function::internals {

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

  //! \brief partial specialization for symmetric tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::stensor<N, real>> {
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::tensor<N, real>> {
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for fourth order tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::t2tot2<N, real>> {
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value *
                                       tfel::math::TensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for fourth order tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::t2tost2<N, real>> {
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value *
                                       tfel::math::StensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for fourth order tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::st2tot2<N, real>> {
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value *
                                       tfel::math::TensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for fourth order tensors
  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::st2tost2<N, real>> {
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value *
                                       tfel::math::StensorDimeToSize<N>::value;
  };

  //! \brief partial specialization for mutable and immutable views
  template <typename T>
  struct CompileTimeSize<tfel::math::View<T>>
      : CompileTimeSize<std::remove_cv_t<T>> {};

  template <typename T>
  struct IsTensor : std::false_type {};

  template <unsigned short N>
  struct IsTensor<tfel::math::fsarray<N, real>> : std::true_type {};

  template <unsigned short N>
  struct IsTensor<tfel::math::tvector<N, real>> : std::true_type {};

  template <unsigned short N, unsigned short M>
  struct IsTensor<tfel::math::tmatrix<N, M, real>> : std::true_type {};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::stensor<N, real>> : std::true_type {
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::tensor<N, real>> : std::true_type {
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::t2tot2<N, real>> : std::true_type {
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::t2tost2<N, real>> : std::true_type {
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::st2tot2<N, real>> : std::true_type {
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::st2tost2<N, real>> : std::true_type {
  };

}  // namespace mgis::function::internals

namespace mgis::function {

  template <typename T>
  concept TensorConcept = internals::IsTensor<T>::value;
  /*!
   * \brief a tensor which is true if the type is equal to real or matches the
   * TensorConcept
   */
  template <typename T>
  concept ScalarOrTensorConcept = std::same_as<T, real> || TensorConcept<T>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORCONCEPT_HXX */
