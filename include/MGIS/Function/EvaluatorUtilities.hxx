/*!
 * \file   MGIS/Function/EvaluatorUtilities.hxx
 * \brief  This file declares a set of utilities for evaluators, such as
 * `evaluator_result`, `number_of_components`, etc..
 * \author Thomas Helfer
 * \date 09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORUTILITIES_HXX
#define LIB_MGIS_FUNCTION_EVALUATORUTILITIES_HXX

#include <span>
#include <array>
#include <concepts>
#include <type_traits>
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function::internals {

  template <bool, EvaluatorConcept EvaluatorType>
  struct EvaluatorResult1 {
    struct InvalidResult;
    using type = InvalidResult;
  };

  template <EvaluatorConcept EvaluatorType>
  struct EvaluatorResult1<true, EvaluatorType> {
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    using type =
        std::invoke_result_t<EvaluatorType, typename Space::element_index_type>;
  };

  template <bool, EvaluatorConcept EvaluatorType>
  struct EvaluatorResult2 {
    struct InvalidResult;
    using type = InvalidResult;
  };

  template <EvaluatorConcept EvaluatorType>
  struct EvaluatorResult2<true, EvaluatorType> {
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    using type = std::invoke_result_t<EvaluatorType,
                                      typename Space::ElementWorkspace,
                                      typename Space::element_index_type>;
  };

  template <bool, EvaluatorConcept EvaluatorType>
  struct EvaluatorResult3 {
    struct InvalidResult;
    using type = InvalidResult;
  };

  template <EvaluatorConcept EvaluatorType>
  struct EvaluatorResult3<true, EvaluatorType> {
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    using type =
        std::invoke_result_t<EvaluatorType,
                             typename Space::cell_index_type,
                             typename Space::quadrature_point_index_type>;
  };

  template <bool, EvaluatorConcept EvaluatorType>
  struct EvaluatorResult4 {
    struct InvalidResult;
    using type = InvalidResult;
  };

  template <EvaluatorConcept EvaluatorType>
  struct EvaluatorResult4<true, EvaluatorType> {
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    using type =
        std::invoke_result_t<EvaluatorType,
                             typename Space::CellWorkspace,
                             typename Space::cell_index_type,
                             typename Space::quadrature_point_index_type>;
  };

  template <EvaluatorConcept EvaluatorType>
  struct EvaluatorResult {
   private:
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    static constexpr bool b1 = requires(EvaluatorType & e) {
      e(std::declval<typename Space::element_index_type>());
    };
    static constexpr bool b2 = requires(EvaluatorType & e) {
      e(std::declval<typename Space::ElementWorkspace>(),
        std::declval<typename Space::element_index_type>());
    };
    static constexpr bool b3 = requires(EvaluatorType & e) {
      e(std::declval<typename Space::CellWorkspace>(),
        std::declval<typename Space::cell_index_type>());
    };
    static constexpr bool b4 = requires(EvaluatorType & e) {
      e(std::declval<typename Space::CellWorkspace>(),
        std::declval<typename Space::cell_index_type>(),
        std::declval<typename Space::quadrature_point_index_type>());
    };
    using ResultType1 = typename EvaluatorResult1<b1, EvaluatorType>::type;
    using ResultType2 = typename EvaluatorResult2<b2, EvaluatorType>::type;
    using ResultType3 = typename EvaluatorResult3<b3, EvaluatorType>::type;
    using ResultType4 = typename EvaluatorResult4<b4, EvaluatorType>::type;
    static_assert(b1 || b2 || b3 || b4);

   public:
    // result of the meta_function
    using type = std::conditional_t<
        b1,
        ResultType1,
        std::conditional_t<b2,
                           ResultType2,
                           std::conditional_t<b3, ResultType3, ResultType4>>>;

   private:
    static_assert(b1 ? std::same_as<ResultType1, type> : true);
    static_assert(b2 ? std::same_as<ResultType2, type> : true);
    static_assert(b3 ? std::same_as<ResultType3, type> : true);
    static_assert(b4 ? std::same_as<ResultType4, type> : true);
  };

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
   * \brief type of the result of an evaluator
   */
  template <EvaluatorConcept EvaluatorType>
  using evaluator_result =
      typename internals::EvaluatorResult<EvaluatorType>::type;

  /*!
   * \brief number of components of a type when known at compile-time,
   * dynamic_extent otherwise
   */
  template <EvaluatorConcept EvaluatorType>
  inline constexpr size_type number_of_components =
      internals::CompileTimeSize<evaluator_result<EvaluatorType>>::value;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORUTILITIES_HXX */
