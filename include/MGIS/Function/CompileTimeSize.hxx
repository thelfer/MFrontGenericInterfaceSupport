/*!
 * \file   MGIS/Function/CompileTimeSize.hxx
 * \brief File which defines the helpers used to determine the number of
 * components of a type (scalar value, span, array, tensor, etc...) when it is
 * known at compile time.
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

#ifndef LIB_MGIS_FUNCTION_COMPILETIMESIZE_HXX
#define LIB_MGIS_FUNCTION_COMPILETIMESIZE_HXX

#include <span>
#include <array>
#include "MGIS/Config.hxx"

namespace mgis::function::internals {

  /*!
   * \brief traits giving the number of components of a type
   *
   * This traits is specialized for scalar values, spans and arrays and, in the
   * `MGIS/Function/TFEL` subdirectory, for TFEL's tensorial types. The primary
   * template applies to the other types, such as `std::vector<real>`, whose
   * number of components is only known at run time.
   *
   * \tparam T: type
   */
  template <typename T>
  struct CompileTimeSize {
    //! \brief default value: the number of components is only known at run time
    static constexpr size_type value = dynamic_extent;
  };

  //! \brief specialization for a scalar value, which has one component
  template <>
  struct CompileTimeSize<real> {
    //! \brief number of components
    static constexpr size_type value = 1;
  };

  /*!
   * \brief specialization for a fixed-sized span
   * \tparam N: number of values of the span
   */
  template <std::size_t N>
  struct CompileTimeSize<std::span<real, N>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  /*!
   * \brief specialization for a constant fixed-sized span
   * \tparam N: number of values of the span
   */
  template <std::size_t N>
  struct CompileTimeSize<std::span<const real, N>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  //! \brief specialization for a dynamic-sized span
  template <>
  struct CompileTimeSize<std::span<real>> {
    //! \brief the number of components is only known at run time
    static constexpr size_type value = dynamic_extent;
  };

  //! \brief specialization for a constant dynamic-sized span
  template <>
  struct CompileTimeSize<std::span<const real>> {
    //! \brief the number of components is only known at run time
    static constexpr size_type value = dynamic_extent;
  };

  /*!
   * \brief specialization for a fixed-sized array
   * \tparam N: number of values of the array
   */
  template <std::size_t N>
  struct CompileTimeSize<std::array<real, N>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  /*!
   * \brief specialization for a constant fixed-sized array
   * \tparam N: number of values of the array
   */
  template <std::size_t N>
  struct CompileTimeSize<std::array<const real, N>> {
    //! \brief number of components
    static constexpr size_type value = N;
  };

  /*!
   * \brief traits giving the number of components of a function or an evaluator
   * when it is known at compile-time, dynamic_extent otherwise
   *
   * \tparam FunctionOrEvaluatorType: function or evaluator type
   * \note this class is specialized for the evaluators in
   * `MGIS/Function/EvaluatorConcept.hxx` and for the functions in
   * `MGIS/Function/FunctionConcept.hxx`
   */
  template <typename FunctionOrEvaluatorType>
  struct NumberOfComponents;

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief number of components (size) of a type if known at compile time,
   * dynamic_extent otherwise
   *
   * \tparam T: type
   */
  template <typename T>
  inline constexpr size_type compile_time_size =
      internals::CompileTimeSize<T>::value;

  /*!
   * \brief number of components of a function or an evaluator when known at
   * compile-time, dynamic_extent otherwise
   *
   * \tparam FunctionOrEvaluatorType: function or evaluator type
   */
  template <typename FunctionOrEvaluatorType>
  inline constexpr size_type number_of_components =
      internals::NumberOfComponents<FunctionOrEvaluatorType>::value;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_COMPILETIMESIZE_HXX */
