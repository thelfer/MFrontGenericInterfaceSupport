/*!
 * \file   MGIS/Function/CompileTimeSize.hxx
 * \brief
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

namespace mgis::function::internals {

  /*!
   * \return the number of components (size) of a type
   * \tparam T: type
   */
  template <typename T>
  struct CompileTimeSize {
    //! \brief default value
    static constexpr size_type value = dynamic_extent;
  };

  //! \brief partial specialization for floating point number
  template <>
  struct CompileTimeSize<real> {
    static constexpr size_type value = 1;
  };

  //! \brief partial specialization for fixed-sized span
  template <std::size_t N>
  struct CompileTimeSize<std::span<real, N>> {
    static constexpr size_type value = N;
  };

  //! \brief partial specialization for constant fixed-sized span
  template <std::size_t N>
  struct CompileTimeSize<std::span<const real, N>> {
    static constexpr size_type value = N;
  };

  //! \brief partial specialization for dynamic-sized span
  template <>
  struct CompileTimeSize<std::span<real>> {
    static constexpr size_type value = dynamic_extent;
  };

  //! \brief partial specialization for constant dynamic-sized span
  template <>
  struct CompileTimeSize<std::span<const real>> {
    static constexpr size_type value = dynamic_extent;
  };

  //! \brief partial specialization for fixed-sized array
  template <std::size_t N>
  struct CompileTimeSize<std::array<real, N>> {
    static constexpr size_type value = N;
  };

  //! \brief partial specialization for constant fixed-sized array
  template <std::size_t N>
  struct CompileTimeSize<std::array<const real, N>> {
    static constexpr size_type value = N;
  };

  /*!
   * \brief number of components of a type when known at compile-time,
   * dynamic_extent otherwise
   *
   * This class is specialized later for evaluators and functions
   */
  template <typename FunctionOrEvaluatorType>
  struct NumberOfComponents;

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \return the number of components (size) of a type if known, dynamic_extent
   * otherwise
   * \tparam T: type
   */
  template <typename T>
  inline constexpr size_type compile_time_size =
      internals::CompileTimeSize<T>::value;

  /*!
   * \brief number of components of a type when known at compile-time,
   * dynamic_extent otherwise
   */
  template <typename FunctionOrEvaluatorType>
  inline constexpr size_type number_of_components =
      internals::NumberOfComponents<FunctionOrEvaluatorType>::value;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_COMPILETIMESIZE_HXX */
