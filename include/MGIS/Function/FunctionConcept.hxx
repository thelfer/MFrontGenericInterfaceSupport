/*!
 * \file   MGIS/Function/FunctionConcept.hxx
 * \brief File which defines the concept a function must satisfy and the
 * helpers used to query the space, the type of the values and the view of a
 * function.
 * \author Thomas Helfer
 * \date   14/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_HXX
#define LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_HXX

#include <span>
#include <array>
#include <utility>
#include <concepts>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function {

  namespace internals {

    /*!
     * \brief traits that must be specialized to allow functions to return this
     * type.
     *
     * \tparam T: type for which the traits are queried
     */
    template <typename T>
    struct FunctionResultTypeTraits {
      //! \brief the traits are not specialized for the queried type
      static constexpr auto is_specialized = false;
    };

    //! \brief partial specialization for reference to a scalar value
    template <>
    struct FunctionResultTypeTraits<real&> {
      //! \brief the traits are specialized for this type
      static constexpr auto is_specialized = true;
    };

    /*!
     * \brief partial specialization for a mutable span
     *
     * \tparam N: number of components of the values returned by the function,
     * which may be `std::dynamic_extent`
     */
    template <std::size_t N>
    struct FunctionResultTypeTraits<std::span<real, N>> {
      //! \brief the traits are specialized for this type
      static constexpr auto is_specialized = true;
    };

    /*!
     * \brief partial specialization for a reference to a mutable array
     *
     * \tparam N: number of components of the values returned by the function
     */
    template <std::size_t N>
    struct FunctionResultTypeTraits<std::array<real, N>&> {
      //! \brief the traits are specialized for this type
      static constexpr auto is_specialized = true;
    };

    /*!
     * \brief concept defining the types a function may return
     *
     * A function may return a mutable scalar value, a mutable `std::span` or a
     * mutable `std::array`, i.e. the types for which the
     * `FunctionResultTypeTraits` traits is specialized.
     *
     * \tparam T: type to be tested
     */
    template <typename T>
    concept mutable_return_value = FunctionResultTypeTraits<T>::is_specialized;

    /*!
     * \brief helper class used to determine the type returned by a function
     * called with an element index, which is an invalid type if the function
     * cannot be called this way
     *
     * \tparam b: false if the function cannot be called with an element index
     * \tparam FunctionType: function type
     */
    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation1 {
      //! \brief type returned when the function cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the function
      using result_type = InvalidResult;
      //! \brief type of the value returned by a constant function
      using const_result_type = InvalidResult;
    };

    /*!
     * \brief specialization used when the function can be called this way
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    struct FunctionResultQueryImplementation1<true, FunctionType> {
      //! \brief space on which the function is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      //! \brief type of the value returned by the function
      using result_type =
          std::invoke_result_t<FunctionType, element_index<Space>>;
      //! \brief type of the value returned by a constant function
      using const_result_type =
          std::invoke_result_t<const FunctionType, element_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by a function
     * called with an element workspace and an element index, which is an
     * invalid type if the function cannot be called this way
     *
     * \tparam b: false if the function cannot be called this way
     * \tparam FunctionType: function type
     */
    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation2 {
      //! \brief type returned when the function cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the function
      using result_type = InvalidResult;
      //! \brief type of the value returned by a constant function
      using const_result_type = InvalidResult;
    };

    /*!
     * \brief specialization used when the function can be called this way
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    struct FunctionResultQueryImplementation2<true, FunctionType> {
      //! \brief space on which the function is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      //! \brief type of the value returned by the function
      using result_type = std::invoke_result_t<FunctionType,
                                               element_workspace<Space>,
                                               element_index<Space>>;
      //! \brief type of the value returned by a constant function
      using const_result_type = std::invoke_result_t<const FunctionType,
                                                     element_workspace<Space>,
                                                     element_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by a function
     * called with a cell index and a quadrature point index, which is an
     * invalid type if the function cannot be called this way
     *
     * \tparam b: false if the function cannot be called this way
     * \tparam FunctionType: function type
     */
    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation3 {
      //! \brief type returned when the function cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the function
      using result_type = InvalidResult;
      //! \brief type of the value returned by a constant function
      using const_result_type = InvalidResult;
    };

    /*!
     * \brief specialization used when the function can be called this way
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    struct FunctionResultQueryImplementation3<true, FunctionType> {
      //! \brief space on which the function is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      //! \brief type of the value returned by the function
      using result_type = std::invoke_result_t<FunctionType,
                                               cell_index<Space>,
                                               quadrature_point_index<Space>>;
      //! \brief type of the value returned by a constant function
      using const_result_type =
          std::invoke_result_t<const FunctionType,
                               cell_index<Space>,
                               quadrature_point_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by a function
     * called with a cell workspace, a cell index and a quadrature point index,
     * which is an invalid type if the function cannot be called this way
     *
     * \tparam b: false if the function cannot be called this way
     * \tparam FunctionType: function type
     */
    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation4 {
      //! \brief type returned when the function cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the function
      using result_type = InvalidResult;
      //! \brief type of the value returned by a constant function
      using const_result_type = InvalidResult;
    };

    /*!
     * \brief specialization used when the function can be called this way
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    struct FunctionResultQueryImplementation4<true, FunctionType> {
      //! \brief space on which the function is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      //! \brief type of the value returned by the function
      using result_type = std::invoke_result_t<FunctionType,
                                               cell_workspace<Space>,
                                               cell_index<Space>,
                                               quadrature_point_index<Space>>;
      //! \brief type of the value returned by a constant function
      using const_result_type =
          std::invoke_result_t<const FunctionType,
                               cell_workspace<Space>,
                               cell_index<Space>,
                               quadrature_point_index<Space>>;
    };

    /*!
     * \brief helper class used to determine how a function can be called and
     * the type of the values it returns
     *
     * \tparam b: false if the function does not provide a `getSpace` function
     * \tparam FunctionType: function type
     */
    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation {
      //! \brief type used for the queries which failed
      struct InvalidResult;
      //! \brief the function can be called with an element index
      static constexpr bool b1 = false;
      //! \brief the function can be called with an element workspace and an
      /// element index
      static constexpr bool b2 = false;
      //! \brief the function can be called with a cell index and a quadrature
      /// point index
      static constexpr bool b3 = false;
      //! \brief the function can be called with a cell workspace, a cell index
      /// and a quadrature point index
      static constexpr bool b4 = false;
      //! \brief type returned by a function called with an element index
      using result_type1 = InvalidResult;
      //! \brief type returned by a function called with an element workspace
      /// and an element index
      using result_type2 = InvalidResult;
      //! \brief type returned by a function called with a cell index and a
      /// quadrature point index
      using result_type3 = InvalidResult;
      //! \brief type returned by a function called with a cell workspace, a
      /// cell index and a quadrature point index
      using result_type4 = InvalidResult;
      //! \brief type returned by a constant function called with an element
      /// index
      using const_result_1 = InvalidResult;
      //! \brief type returned by a constant function called with an element
      /// workspace and an element index
      using const_result_2 = InvalidResult;
      //! \brief type returned by a constant function called with a cell index
      /// and a quadrature point index
      using const_result_3 = InvalidResult;
      //! \brief type returned by a constant function called with a cell
      /// workspace, a cell index and a quadrature point index
      using const_result_4 = InvalidResult;
      //! \brief type of the values returned by the function
      using result_type = InvalidResult;
      //! \brief type of the values returned by a constant function
      using const_result_type = InvalidResult;
    };

    /*!
     * \brief specialization used when the function provides a `getSpace`
     * function
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    struct FunctionResultQueryImplementation<true, FunctionType> {
      //! \brief space on which the function is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      //! \brief the function can be called with an element index
      static constexpr bool b1 = ((requires(FunctionType & e) {
                                    {
                                      e(std::declval<element_index<Space>>())
                                      } -> mutable_return_value;
                                  }) &&
                                  (ElementSpaceConcept<Space>));
      //! \brief the function can be called with an element workspace and an
      /// element index
      static constexpr bool b2 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<element_workspace<Space>>(),
                 std::declval<element_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (ElementSpaceConcept<Space> && hasElementWorkspace<Space>));
      //! \brief the function can be called with a cell index and a quadrature
      /// point index
      static constexpr bool b3 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<cell_index<Space>>(),
                 std::declval<quadrature_point_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (QuadratureSpaceConcept<Space>));
      //! \brief the function can be called with a cell workspace, a cell index
      /// and a quadrature point index
      static constexpr bool b4 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<cell_workspace<Space>>(),
                 std::declval<cell_index<Space>>(),
                 std::declval<quadrature_point_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (QuadratureSpaceConcept<Space> && hasCellWorkspace<Space>));
      //! \brief type returned by a function called with an element index
      using result_type1 =
          typename FunctionResultQueryImplementation1<b1, FunctionType>::
              result_type;
      //! \brief type returned by a function called with an element workspace
      /// and an element index
      using result_type2 =
          typename FunctionResultQueryImplementation2<b2, FunctionType>::
              result_type;
      //! \brief type returned by a function called with a cell index and a
      /// quadrature point index
      using result_type3 =
          typename FunctionResultQueryImplementation3<b3, FunctionType>::
              result_type;
      //! \brief type returned by a function called with a cell workspace, a
      /// cell index and a quadrature point index
      using result_type4 =
          typename FunctionResultQueryImplementation4<b4, FunctionType>::
              result_type;
      //! \brief type of the values returned by the function, which is the type
      /// returned by the first valid call
      using result_type = std::conditional_t<
          b1,
          result_type1,
          std::conditional_t<
              b2,
              result_type2,
              std::conditional_t<b3, result_type3, result_type4>>>;
      //! \brief type returned by a constant function called with an element
      /// index
      using const_result_type1 =
          typename FunctionResultQueryImplementation1<b1, FunctionType>::
              const_result_type;
      //! \brief type returned by a constant function called with an element
      /// workspace and an element index
      using const_result_type2 =
          typename FunctionResultQueryImplementation2<b2, FunctionType>::
              const_result_type;
      //! \brief type returned by a constant function called with a cell index
      /// and a quadrature point index
      using const_result_type3 =
          typename FunctionResultQueryImplementation3<b3, FunctionType>::
              const_result_type;
      //! \brief type returned by a constant function called with a cell
      /// workspace, a cell index and a quadrature point index
      using const_result_type4 =
          typename FunctionResultQueryImplementation4<b4, FunctionType>::
              const_result_type;
      //! \brief type of the values returned by a constant function, which is
      /// the type returned by the first valid call
      using const_result_type = std::conditional_t<
          b1,
          const_result_type1,
          std::conditional_t<
              b2,
              const_result_type2,
              std::conditional_t<b3, const_result_type3, const_result_type4>>>;
    };

    /*!
     * \brief alias used to determine how a function can be called and the type
     * of the values it returns
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
        struct FunctionResultQuery
        : FunctionResultQueryImplementation < requires(FunctionType& e) {
      getSpace(e);
    }, FunctionType > {};

  }  // namespace internals

  /*!
   * \brief traits defining the functions which are their own view
   *
   * This traits must be specialized for a function type to indicate that the
   * function is lightweight, i.e. that a copy of the function is cheap and can
   * be used to evaluate the function without creating a view.
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  struct LightweightViewTraits : std::false_type {};

  /*!
   * \brief concept defining the functions which are their own view
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept LightweightViewConcept = LightweightViewTraits<FunctionType>::value;

  /*!
   * \brief a concept that must satisfy a function
   *
   * A function is an object which:
   * - is defined on a space, i.e. for which the `getSpace` function is
   *   available,
   * - provides the number of components of the values it returns through the
   *   `getNumberOfComponents` function,
   * - can be called with either an element index, or an element workspace and
   *   an element index, or a cell index and a quadrature point index, or a
   *   cell workspace, a cell index and a quadrature point index. The value
   *   returned is either a scalar value, a `std::span` or a `std::array`,
   * - is either a lightweight view or provides a `view` function returning such
   *   a view.
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept FunctionConcept =
      (SpaceConcept<
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>>)&&  //
      (requires(const FunctionType& e) {
        { getNumberOfComponents(e) } -> std::same_as<mgis::size_type>;
      }) &&
      ((internals::FunctionResultQuery<FunctionType>::b1) ||
       (internals::FunctionResultQuery<FunctionType>::b2) ||
       (internals::FunctionResultQuery<FunctionType>::b3) ||
       (internals::FunctionResultQuery<FunctionType>::b4)) &&
      (internals::FunctionResultQuery<FunctionType>::b1
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::result_type1,
                          typename internals::FunctionResultQuery<
                              FunctionType>::result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b2
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::result_type2,
                          typename internals::FunctionResultQuery<
                              FunctionType>::result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b3
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::result_type3,
                          typename internals::FunctionResultQuery<
                              FunctionType>::result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b4
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::result_type4,
                          typename internals::FunctionResultQuery<
                              FunctionType>::result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b1
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type1,
                          typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b2
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type2,
                          typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b3
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type3,
                          typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type>
           : true) &&
      (internals::FunctionResultQuery<FunctionType>::b4
           ? std::same_as<typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type4,
                          typename internals::FunctionResultQuery<
                              FunctionType>::const_result_type>
           : true) &&
      (internals::FunctionResultTypeTraits<
          typename internals::FunctionResultQuery<FunctionType>::result_type>::
           is_specialized) &&
      ((LightweightViewConcept<FunctionType>) || (requires(FunctionType & f) {
         { view(f) } -> LightweightViewConcept;
       }));

  /*!
   * \brief concept defining the functions which are their own view
   *
   * Such a function can be copied at a negligible cost, which makes the
   * creation of a view useless.
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept LightweightFunctionConcept =
      (FunctionConcept<FunctionType>)&&(LightweightViewConcept<FunctionType>);

  /*!
   * \brief concept defining the functions which are not their own view
   *
   * The value of such a function can only be accessed through a view, which is
   * created by the `view` function.
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept NonLightweightFunctionConcept =
      (FunctionConcept<FunctionType>)&&(!LightweightViewConcept<FunctionType>);

  /*!
   * \brief concept defining the arguments of the algorithms which can be
   * viewed as a mutable function
   *
   * A constant function or a function which is not a lvalue cannot be used as a
   * mutable function argument.
   *
   * \tparam QualifiedFunctionArgumentType: type of the function argument, which
   * may be a qualified type (reference, const, etc...)
   */
  template <typename QualifiedFunctionArgumentType>
  concept ViewableFunctionArgumentConcept =
      (FunctionConcept<std::decay_t<QualifiedFunctionArgumentType>>)&&  //
      (!std::is_const_v<
          std::remove_reference_t<QualifiedFunctionArgumentType>>)&&  //
      ((LightweightViewConcept<std::decay_t<QualifiedFunctionArgumentType>>) ||
       (std::is_lvalue_reference_v<QualifiedFunctionArgumentType>));

  namespace internals {

    /*!
     * \brief helper structure used to determine the type of the view of a
     * function which is not its own view
     *
     * \tparam b: true if the function is its own view
     * \tparam FunctionType: function type
     */
    template <bool, FunctionConcept FunctionType>
    struct ViewResultType {
      //! \brief no view is required
      using type = void;
    };

    /*!
     * \brief specialization used when a view is required
     *
     * \tparam FunctionType: function type
     */
    template <FunctionConcept FunctionType>
    struct ViewResultType<false, FunctionType> {
      //! \brief type of the view of the function
      using type = decltype(view(std::declval<FunctionType&>()));
    };

  }  // namespace internals

  /*!
   * \brief an alias which returns a view corresponding to the given function
   * type
   *
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  using function_view = std::conditional_t<
      LightweightViewConcept<FunctionType>,
      FunctionType,
      typename internals::ViewResultType<LightweightViewConcept<FunctionType>,
                                         FunctionType>::type>;

  /*!
   * \brief an helper function to create a view when required
   *
   * \tparam FunctionType: type of the function, which may be a qualified type
   * (reference, const, etc...)
   * \param[in] f: function
   * \return a view of the function, which is the function itself if the
   * function is a lightweight view
   */
  template <ViewableFunctionArgumentConcept FunctionType>
  constexpr function_view<std::decay_t<FunctionType>> make_view(
      FunctionType&& f) {
    if constexpr (LightweightFunctionConcept<std::decay_t<FunctionType>>) {
      return f;
    } else {
      return view(f);
    }
  }  // end of make_view

  /*!
   * \brief concept defining functions which can be called on an element space,
   * i.e. with an element index or with an element workspace and an element
   * index
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept ElementFunctionConcept = FunctionConcept<FunctionType> &&
      ((internals::FunctionResultQuery<FunctionType>::b1) ||
       (internals::FunctionResultQuery<FunctionType>::b2));

  /*!
   * \brief concept defining functions which can be called on a quadrature
   * space, i.e. with a cell index and a quadrature point index or with a cell
   * workspace, a cell index and a quadrature point index
   *
   * \tparam FunctionType: function type
   */
  template <typename FunctionType>
  concept QuadratureFunctionConcept = FunctionConcept<FunctionType> &&
      ((internals::FunctionResultQuery<FunctionType>::b3) ||
       (internals::FunctionResultQuery<FunctionType>::b4));

  /*!
   * \brief alias giving the type of the space on which a function is defined
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  using function_space =
      std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;

  /*!
   * \brief type of the result of a function, i.e. the type of the values
   * returned by the function, which can be a scalar value, a `std::span` or a
   * `std::array`
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  using function_result =
      typename internals::FunctionResultQuery<FunctionType>::result_type;
  /*!
   * \brief type of the result of a constant function, i.e. the type of the
   * values returned by a constant function, which can be a scalar value, a
   * constant `std::span` or a constant `std::array`
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  using function_const_result =
      typename internals::FunctionResultQuery<FunctionType>::const_result_type;

  /*!
   * \brief true if the result of a function is a mappable type, i.e. a type for
   * which the `getDataPointer` function returns a mutable pointer
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  inline constexpr auto isFunctionResultTypeMappable =
      requires(function_result<FunctionType> rf) {
    { getDataPointer(rf) } -> std::same_as<real*>;
  };

  /*!
   * \brief true if the result of a constant function is a mappable type, i.e. a
   * type for which the `getDataPointer` function returns a constant pointer
   * \tparam FunctionType: function type
   */
  template <FunctionConcept FunctionType>
  inline constexpr auto isFunctionConstResultTypeMappable =
      requires(function_const_result<FunctionType> rf) {
    { getDataPointer(rf) } -> std::same_as<const real*>;
  };

  namespace internals {

    /*!
     * \brief This helper function allows to disambiguate the call to
     * the getSpace function
     *
     * \tparam FunctionType: function type
     * \param[in] f: function
     * \return the space on which the function is defined
     */
    template <FunctionConcept FunctionType>
    constexpr decltype(auto) disambiguateGetSpace(const FunctionType& f)  //
        requires(!EvaluatorConcept<FunctionType>);
    /*!
     * \brief This helper function allows to disambiguate the call to
     * the getNumberOfComponents function
     *
     * \tparam FunctionType: function type
     * \param[in] f: function
     * \return the number of components of the values returned by the function
     */
    template <FunctionConcept FunctionType>
    constexpr mgis::size_type disambiguateGetNumberOfComponents(
        const FunctionType& f)  //
        requires(!EvaluatorConcept<FunctionType>);

    /*!
     * \brief number of components of a function when known at compile-time,
     * dynamic_extent otherwise
     *
     * \tparam FunctionType: function type
     */
    template <typename FunctionType>
    requires(FunctionConcept<FunctionType> &&
             (!EvaluatorConcept<FunctionType>))  //
        struct NumberOfComponents<FunctionType> {
      //! \brief number of components of the values returned by the function
      static constexpr size_type value = internals::CompileTimeSize<
          std::decay_t<function_result<FunctionType>>>::value;
    };

    /*!
     * \brief checks if the number of components of a function is compatible
     * with the number of components expected by a modifier
     *
     * The check is always satisfied if the number of components of the function
     * is only known at run time.
     *
     * \tparam QualifiedFunctionArgumentType: type of the function argument,
     * which may be a qualified type (reference, const, etc...)
     * \tparam N: number of components expected by the modifier
     * \return true if the number of components is compatible
     */
    template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType,
              size_type N>
    [[nodiscard]] constexpr bool checkNumberOfComponentsCompatibility()
        MGIS_NOEXCEPT {
      constexpr auto Nc =
          number_of_components<std::decay_t<QualifiedFunctionArgumentType>>;
      if constexpr (Nc == dynamic_extent) {
        return true;
      } else {
        return N == Nc;
      }
    }  // end of checkNumberOfComponentsCompatibility

  }  // namespace internals

#ifndef _MSC_VER
  /*!
   * \brief assign an evaluator to a mutable function
   *
   * \tparam EvaluatorType: evaluator type
   * \tparam FunctionType: type of the function, which may be a qualified type
   * (reference, const, etc...)
   * \param[in] e: evaluator
   * \param[in,out] f: function
   * \return true if the assignment succeeded
   */
  template <EvaluatorConcept EvaluatorType,
            ViewableFunctionArgumentConcept FunctionType>
  [[nodiscard]] bool operator|(EvaluatorType e, FunctionType&& f) requires(
      std::same_as<evaluator_space<EvaluatorType>,
                   function_space<std::decay_t<FunctionType>>>);
#endif

}  // end of namespace mgis::function

#include "MGIS/Function/FunctionConcept.ixx"
#include "MGIS/Function/FixedSizeView.hxx"
#include "MGIS/Function/FixedSizeModifier.hxx"

namespace mgis::function::internals {

  /*!
   * \brief structure used to modify a function or an evaluator so that it
   * returns a fixed number of components
   *
   * \tparam N: number of components
   */
  template <size_type N>
  requires((N > 0) && (N != dynamic_extent))  //
      struct fixed_size_modifier {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;

    /*!
     * \brief create a view returning a fixed number of components
     *
     * \tparam QualifiedFunctionArgumentType: type of the function, which may be
     * a qualified type (reference, const, etc...)
     * \param[in] f: function the view is built from
     *
     * \return a view of the function returning the first `N` components of the
     * values of the function
     */
    template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType>
    constexpr auto operator()(QualifiedFunctionArgumentType&& f) const
        requires(internals::checkNumberOfComponentsCompatibility<
                 QualifiedFunctionArgumentType,
                 N>()) {
      return FixedSizeView<std::decay_t<QualifiedFunctionArgumentType>, N>(
          std::forward<QualifiedFunctionArgumentType>(f));
    }
    /*!
     * \brief create a new evaluator returning a fixed number of components
     *
     * \tparam EvaluatorType: evaluator type
     * \param[in] e: evaluator the new evaluator is built from
     *
     * \return an evaluator returning the first `N` components of the values of
     * the given evaluator
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType& e) const requires(
        checkEvaluatorNumberOfComponentsCompatibility<EvaluatorType, N>()) {
      return FixedSizeModifier<EvaluatorType, N>(e);
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief create a view returning a fixed number of components
   *
   * \tparam QualifiedFunctionArgumentType: type of the function, which may be a
   * qualified type (reference, const, etc...)
   * \tparam N: number of components
   * \param[in] f: function
   * \param[in] m: modifier
   * \return a view of the function returning the first `N` components of the
   * values of the function
   */
  template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType,
            size_type N>
  constexpr auto operator|(QualifiedFunctionArgumentType&& f,
                           const internals::fixed_size_modifier<N>& m)  //
      requires(internals::checkNumberOfComponentsCompatibility<
               QualifiedFunctionArgumentType,
               N>()) {
    return m(std::forward<QualifiedFunctionArgumentType>(f));
  }

  /*!
   * \brief modifier used to return a fixed number of components of a function
   * or an evaluator
   *
   * \tparam N: number of components
   */
  template <size_type N>
  requires((N > 0) && (N != dynamic_extent))  //
      inline constexpr auto as_array = internals::fixed_size_modifier<N>{};

  //! \brief modifier used to return the first component of a function or an
  /// evaluator
  inline constexpr auto as_scalar = internals::fixed_size_modifier<1>{};

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_HXX */
