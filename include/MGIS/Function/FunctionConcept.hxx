/*!
 * \file   MGIS/Function/FunctionConcept.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2025
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
     */
    template <typename T>
    struct FunctionResultTypeTraits {
      static constexpr auto is_specialized = false;
    };

    //! \brief partial specialization for reference to a scalar value
    template <>
    struct FunctionResultTypeTraits<real&> {
      static constexpr auto is_specialized = true;
    };

    //! \brief partial specialization for reference to a mutable span
    template <std::size_t N>
    struct FunctionResultTypeTraits<std::span<real, N>> {
      static constexpr auto is_specialized = true;
    };

    //! \brief partial specialization for reference to a std::array
    template <std::size_t N>
    struct FunctionResultTypeTraits<std::array<real, N>&> {
      static constexpr auto is_specialized = true;
    };

    template <typename T>
    concept mutable_return_value = FunctionResultTypeTraits<T>::is_specialized;

    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation1 {
      struct InvalidResult;
      using result_type = InvalidResult;
      using const_result_type = InvalidResult;
    };

    template <typename FunctionType>
    struct FunctionResultQueryImplementation1<true, FunctionType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      using result_type =
          std::invoke_result_t<FunctionType, element_index<Space>>;
      using const_result_type =
          std::invoke_result_t<const FunctionType, element_index<Space>>;
    };

    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation2 {
      struct InvalidResult;
      using result_type = InvalidResult;
      using const_result_type = InvalidResult;
    };

    template <typename FunctionType>
    struct FunctionResultQueryImplementation2<true, FunctionType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      using result_type = std::invoke_result_t<FunctionType,
                                               element_workspace<Space>,
                                               element_index<Space>>;
      using const_result_type = std::invoke_result_t<const FunctionType,
                                                     element_workspace<Space>,
                                                     element_index<Space>>;
    };

    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation3 {
      struct InvalidResult;
      using result_type = InvalidResult;
      using const_result_type = InvalidResult;
    };

    template <typename FunctionType>
    struct FunctionResultQueryImplementation3<true, FunctionType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      using const_result_type =
          std::invoke_result_t<const FunctionType,
                               cell_index<Space>,
                               quadrature_point_index<Space>>;
    };

    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation4 {
      struct InvalidResult;
      using result_type = InvalidResult;
      using const_result_type = InvalidResult;
    };

    template <typename FunctionType>
    struct FunctionResultQueryImplementation4<true, FunctionType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      using result_type = std::invoke_result_t<FunctionType,
                                               cell_workspace<Space>,
                                               cell_index<Space>,
                                               quadrature_point_index<Space>>;
      using const_result_type =
          std::invoke_result_t<const FunctionType,
                               cell_workspace<Space>,
                               cell_index<Space>,
                               quadrature_point_index<Space>>;
    };

    template <bool, typename FunctionType>
    struct FunctionResultQueryImplementation {
      struct InvalidResult;
      static constexpr bool b1 = false;
      static constexpr bool b2 = false;
      static constexpr bool b3 = false;
      static constexpr bool b4 = false;
      using result_type1 = InvalidResult;
      using result_type2 = InvalidResult;
      using result_type3 = InvalidResult;
      using result_type4 = InvalidResult;
      using const_result_1 = InvalidResult;
      using const_result_2 = InvalidResult;
      using const_result_3 = InvalidResult;
      using const_result_4 = InvalidResult;
      using result_type = InvalidResult;
      using const_result_type = InvalidResult;
    };

    template <typename FunctionType>
    struct FunctionResultQueryImplementation<true, FunctionType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;
      static constexpr bool b1 = ((requires(FunctionType & e) {
                                    {
                                      e(std::declval<element_index<Space>>())
                                      } -> mutable_return_value;
                                  }) &&
                                  (ElementSpaceConcept<Space>));
      static constexpr bool b2 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<element_workspace<Space>>(),
                 std::declval<element_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (ElementSpaceConcept<Space> && hasElementWorkspace<Space>));
      static constexpr bool b3 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<cell_index<Space>>(),
                 std::declval<quadrature_point_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (QuadratureSpaceConcept<Space>));
      static constexpr bool b4 =
          ((requires(FunctionType & e) {
             {
               e(std::declval<cell_workspace<Space>>(),
                 std::declval<cell_index<Space>>(),
                 std::declval<quadrature_point_index<Space>>())
               } -> mutable_return_value;
           }) &&
           (QuadratureSpaceConcept<Space> && hasCellWorkspace<Space>));
      using result_type1 =
          typename FunctionResultQueryImplementation1<b1, FunctionType>::
              result_type;
      using result_type2 =
          typename FunctionResultQueryImplementation2<b2, FunctionType>::
              result_type;
      using result_type3 =
          typename FunctionResultQueryImplementation3<b3, FunctionType>::
              result_type;
      using result_type4 =
          typename FunctionResultQueryImplementation4<b4, FunctionType>::
              result_type;
      using result_type = std::conditional_t<
          b1,
          result_type1,
          std::conditional_t<
              b2,
              result_type2,
              std::conditional_t<b3, result_type3, result_type4>>>;
      using const_result_type1 =
          typename FunctionResultQueryImplementation1<b1, FunctionType>::
              const_result_type;
      using const_result_type2 =
          typename FunctionResultQueryImplementation2<b2, FunctionType>::
              const_result_type;
      using const_result_type3 =
          typename FunctionResultQueryImplementation3<b3, FunctionType>::
              const_result_type;
      using const_result_type4 =
          typename FunctionResultQueryImplementation4<b4, FunctionType>::
              const_result_type;
      using const_result_type = std::conditional_t<
          b1,
          const_result_type1,
          std::conditional_t<
              b2,
              const_result_type2,
              std::conditional_t<b3, const_result_type3, const_result_type4>>>;
    };

    template <typename FunctionType>
        struct FunctionResultQuery
        : FunctionResultQueryImplementation < requires(FunctionType& e) {
      getSpace(e);
    }, FunctionType > {};

  }  // namespace internals

  template <typename FunctionType>
  struct LightweightViewTraits : std::false_type {};

  template <typename FunctionType>
  concept LightweightViewConcept = LightweightViewTraits<FunctionType>::value;

  //! \brief a concept that must satisfy a function
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

  namespace internals {

    template <bool, FunctionConcept FunctionType>
    struct ViewResultType {
      using type = void;
    };

    template <FunctionConcept FunctionType>
    struct ViewResultType<false, FunctionType> {
      using type = decltype(view(std::declval<FunctionType&>()));
    };

  }  // namespace internals

  /*!
   * \brief an alias which returns a view corresponding to the given function
   * type
   */
  template <FunctionConcept FunctionType>
  using function_view = std::conditional_t<
      LightweightViewConcept<FunctionType>,
      FunctionType,
      typename internals::ViewResultType<LightweightViewConcept<FunctionType>,
                                         FunctionType>::type>;

  //! \brief an helper function to create a view when required
  template <FunctionConcept FunctionType>
  constexpr function_view<FunctionType> make_view(FunctionType& f) {
    if constexpr (LightweightViewConcept<FunctionType>) {
      return f;
    } else {
      return view(f);
    }
  }  // end of make_view

  //! \brief concept defining evaluators working on an element space
  template <typename FunctionType>
  concept ElementFunctionConcept = FunctionConcept<FunctionType> &&
      ((internals::FunctionResultQuery<FunctionType>::b1) ||
       (internals::FunctionResultQuery<FunctionType>::b2));

  //! \brief concept defining evaluators working on a quadrature space
  template <typename FunctionType>
  concept QuadratureFunctionConcept = FunctionConcept<FunctionType> &&
      ((internals::FunctionResultQuery<FunctionType>::b3) ||
       (internals::FunctionResultQuery<FunctionType>::b4));

  //
  template <FunctionConcept FunctionType>
  using function_space =
      std::decay_t<decltype(getSpace(std::declval<FunctionType>()))>;

  /*!
   * \brief type of the result of an function
   */
  template <FunctionConcept FunctionType>
  using function_result =
      typename internals::FunctionResultQuery<FunctionType>::result_type;
  /*!
   * \brief type of the result of an function
   */
  template <FunctionConcept FunctionType>
  using function_const_result =
      typename internals::FunctionResultQuery<FunctionType>::const_result_type;

  template <FunctionConcept FunctionType>
  inline constexpr auto isFunctionResultTypeMappable =
      requires(function_result<FunctionType> rf) {
    { rf.data() } -> std::same_as<real*>;
  };

  template <FunctionConcept FunctionType>
  inline constexpr auto isFunctionConstResultTypeMappable =
      requires(function_const_result<FunctionType> rf) {
    { rf.data() } -> std::same_as<const real*>;
  };

  namespace internals {

    /*!
     * This helper function allows to disambiguate the call to
     * the getSpace function
     */
    template <FunctionConcept FunctionType>
    constexpr decltype(auto) disambiguateGetSpace(const FunctionType&)  //
        requires(!EvaluatorConcept<FunctionType>);
    /*!
     * This helper function allows to disambiguate the call to
     * the getNumberOfComponents function
     */
    template <FunctionConcept FunctionType>
    constexpr mgis::size_type disambiguateGetNumberOfComponents(
        const FunctionType&)  //
        requires(!EvaluatorConcept<FunctionType>);

    /*!
     * \brief number of components of a type when known at compile-time,
     * dynamic_extent otherwise
     */
    template <typename FunctionType>
    requires(FunctionConcept<FunctionType> &&
             (!EvaluatorConcept<FunctionType>))  //
        struct NumberOfComponents<FunctionType> {
      static constexpr size_type value = internals::CompileTimeSize<
          std::decay_t<function_result<FunctionType>>>::value;
    };

  }  // namespace internals

#ifndef _MSC_VER
  /*!
   * \brief assign an evaluator to a mutable function view
   * \param[in] ctx: execution context
   * \param[in] e: evaluator
   * \param[in] f: function
   */
  template <EvaluatorConcept EvaluatorType, FunctionConcept FunctionType>
  [[nodiscard]] bool operator|(EvaluatorType, FunctionType&) requires(
      std::same_as<evaluator_space<EvaluatorType>,
                   function_space<FunctionType>>);
#endif

}  // end of namespace mgis::function

#include "MGIS/Function/FunctionConcept.ixx"
#include "MGIS/Function/FixedSizeView.hxx"
#include "MGIS/Function/FixedSizeModifier.hxx"

namespace mgis::function::internals {

  template <size_type N>
  requires((N > 0) && (N != dynamic_extent))  //
      struct fixed_size_modifier {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;

    template <FunctionConcept FunctionType>
    constexpr auto operator()(FunctionType& f) const
        requires(number_of_components<FunctionType> == dynamic_extent
                     ? true
                     : N == number_of_components<FunctionType>) {
      return FixedSizeView<FunctionType, N>(f);
    }
    /*!
     * \brief create a new modifier
     * \param[in] e: evaluator type
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType& e) const
        requires(number_of_components<EvaluatorType> == dynamic_extent
                     ? true
                     : N == number_of_components<EvaluatorType>) {
      return FixedSizeModifier<EvaluatorType, N>(e);
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <FunctionConcept FunctionType, size_type N>
  constexpr auto operator|(FunctionType& f,
                           const internals::fixed_size_modifier<N>& m)  //
      requires(number_of_components<FunctionType> == dynamic_extent
                   ? true
                   : N == number_of_components<FunctionType>) {
    return m(f);
  }

  template <size_type N>
  requires((N > 0) && (N != dynamic_extent))  //
      inline constexpr auto as_array = internals::fixed_size_modifier<N>{};

  inline constexpr auto as_scalar = internals::fixed_size_modifier<1>{};

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_HXX */
