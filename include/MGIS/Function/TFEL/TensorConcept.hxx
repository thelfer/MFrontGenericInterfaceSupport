/*!
 * \file   MGIS/Function/TFEL/TensorConcept.hxx
 * \brief  File which defines the concepts describing TFEL's tensorial values
 * and the specialisations of the `CompileTimeSize` and `IsTensor` helpers for
 * TFEL's fixed size tensorial types.
 * \author Thomas Helfer
 * \date   10/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use tensor evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORCONCEPT_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORCONCEPT_HXX

#include <concepts>
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
#include "MGIS/Function/TFEL/General.hxx"

namespace mgis::function {

  /*!
   * \brief concept defining the scalar types which may be used as value type
   * of a tensorial object
   *
   * A type matches this concept if it matches the
   * `tfel::math::MutableScalarConcept` concept of TFEL and if its base
   * numeric type, as returned by `tfel::math::base_type`, is `real`.
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept TensorValueConcept = ((::tfel::math::MutableScalarConcept<T>)&&(
      std::same_as<::tfel::math::base_type<T>, real>));

}  // end of namespace mgis::function

namespace mgis::function::internals {

  /*!
   * \brief partial specialization for finite size arrays
   *
   * \tparam N: number of values stored in the array
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  struct CompileTimeSize<tfel::math::fsarray<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  /*!
   * \brief partial specialization for tiny vectors
   *
   * \tparam N: number of components of the vector
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  struct CompileTimeSize<tfel::math::tvector<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  /*!
   * \brief partial specialization for tiny matrices
   *
   * \tparam N: number of rows of the matrix
   * \tparam M: number of columns of the matrix
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, unsigned short M, TensorValueConcept ValueType>
  struct CompileTimeSize<tfel::math::tmatrix<N, M, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = N * M;
  };

  /*!
   * \brief partial specialization for symmetric tensors of dimension `N`,
   * where `N` is equal to 1, 2 or 3
   *
   * \tparam N: dimension of the symmetric tensor
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::stensor<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for tensors of dimension `N`, where `N` is
   * equal to 1, 2 or 3
   *
   * \tparam N: dimension of the tensor
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::tensor<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for fourth order tensors describing a
   * linear mapping from tensors to tensors, of dimension `N`, where `N` is
   * equal to 1, 2 or 3
   *
   * \tparam N: dimension of the tensors mapped
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::t2tot2<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value *
                                       tfel::math::TensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for fourth order tensors describing a
   * linear mapping from tensors to symmetric tensors, of dimension `N`, where
   * `N` is equal to 1, 2 or 3
   *
   * \tparam N: dimension of the tensors mapped
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::t2tost2<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::TensorDimeToSize<N>::value *
                                       tfel::math::StensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for fourth order tensors describing a
   * linear mapping from symmetric tensors to tensors, of dimension `N`, where
   * `N` is equal to 1, 2 or 3
   *
   * \tparam N: dimension of the tensors mapped
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::st2tot2<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value *
                                       tfel::math::TensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for fourth order tensors describing a
   * linear mapping from symmetric tensors to symmetric tensors, of dimension
   * `N`, where `N` is equal to 1, 2 or 3
   *
   * \tparam N: dimension of the symmetric tensors mapped
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CompileTimeSize<tfel::math::st2tost2<N, ValueType>> {
    //! \brief number of components
    static constexpr size_type value = tfel::math::StensorDimeToSize<N>::value *
                                       tfel::math::StensorDimeToSize<N>::value;
  };

  /*!
   * \brief partial specialization for mutable and immutable views
   *
   * The number of components of a view is that of the viewed type, constant
   * qualifiers removed.
   *
   * \tparam T: type of the viewed object, which may be constant
   */
  template <typename T>
  struct CompileTimeSize<tfel::math::View<T>>
      : CompileTimeSize<std::remove_cv_t<T>> {};

  /*!
   * \brief traits stating if a type may be used as a tensorial object
   *
   * The primary template states that the type is not a tensorial object. It is
   * partially specialized below for TFEL's fixed size tensorial types, which
   * can all be mapped to a sequence of values stored contiguously in memory.
   *
   * \tparam T: type considered
   */
  template <typename T>
  struct IsTensor : std::false_type {};

  //! \brief specialization for finite size arrays
  template <unsigned short N, TensorValueConcept ValueType>
  struct IsTensor<tfel::math::fsarray<N, ValueType>> : std::true_type {};

  //! \brief specialization for tiny vectors
  template <unsigned short N, TensorValueConcept ValueType>
  struct IsTensor<tfel::math::tvector<N, ValueType>> : std::true_type {};

  //! \brief specialization for tiny matrices
  template <unsigned short N, unsigned short M, TensorValueConcept ValueType>
  struct IsTensor<tfel::math::tmatrix<N, M, ValueType>> : std::true_type {};

  /*!
   * \brief specialization for symmetric tensors whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::stensor<N, ValueType>> : std::true_type {
  };

  /*!
   * \brief specialization for tensors whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::tensor<N, ValueType>> : std::true_type {
  };

  /*!
   * \brief specialization for fourth order tensors mapping tensors to
   * tensors, whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::t2tot2<N, ValueType>> : std::true_type {
  };

  /*!
   * \brief specialization for fourth order tensors mapping tensors to
   * symmetric tensors, whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::t2tost2<N, ValueType>> : std::true_type {
  };

  /*!
   * \brief specialization for fourth order tensors mapping symmetric tensors
   * to tensors, whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::st2tot2<N, ValueType>> : std::true_type {
  };

  /*!
   * \brief specialization for fourth order tensors mapping symmetric tensors
   * to symmetric tensors, whose dimension is 1, 2 or 3
   */
  template <unsigned short N, TensorValueConcept ValueType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct IsTensor<tfel::math::st2tost2<N, ValueType>> : std::true_type {
  };

}  // namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief concept defining the types which may be used as tensorial objects
   *
   * A type matches this concept if the `internals::IsTensor` traits is
   * specialized for it, i.e. if it is one of TFEL's fixed size tensorial
   * types: finite size array, tiny vector, tiny matrix, symmetric tensor,
   * tensor or fourth order tensor.
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept TensorConcept = internals::IsTensor<T>::value;

  /*!
   * \brief concept which is satisfied if a type is equal to `real` or matches
   * the `TensorConcept` concept
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept ScalarOrTensorConcept = TensorValueConcept<T> || TensorConcept<T>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORCONCEPT_HXX */
