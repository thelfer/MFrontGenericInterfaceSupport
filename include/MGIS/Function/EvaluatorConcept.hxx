/*!
 * \file   MGIS/Function/EvaluatorConcept.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX
#define LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX

#include <concepts>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis {

  // forward declaration
  struct AbstractErrorHandler;

}  // end of namespace mgis

namespace mgis::function {

  namespace internals {

    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation1 {
      struct InvalidResult;
      using type = InvalidResult;
    };

    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation1<true, EvaluatorType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      using type =
          std::invoke_result_t<const EvaluatorType, element_index<Space>>;
    };

    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation2 {
      struct InvalidResult;
      using type = InvalidResult;
    };

    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation2<true, EvaluatorType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      using type = std::invoke_result_t<const EvaluatorType,
                                        element_workspace<Space>,
                                        element_index<Space>>;
    };

    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation3 {
      struct InvalidResult;
      using type = InvalidResult;
    };

    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation3<true, EvaluatorType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      using type = std::invoke_result_t<const EvaluatorType,
                                        cell_index<Space>,
                                        quadrature_point_index<Space>>;
    };

    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation4 {
      struct InvalidResult;
      using type = InvalidResult;
    };

    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation4<true, EvaluatorType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      using type = std::invoke_result_t<const EvaluatorType,
                                        cell_workspace<Space>,
                                        cell_index<Space>,
                                        quadrature_point_index<Space>>;
    };

    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation {
      struct InvalidResult;
      static constexpr bool b1 = false;
      static constexpr bool b2 = false;
      static constexpr bool b3 = false;
      static constexpr bool b4 = false;
      using ResultType1 = InvalidResult;
      using ResultType2 = InvalidResult;
      using ResultType3 = InvalidResult;
      using ResultType4 = InvalidResult;
    };

    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation<true, EvaluatorType> {
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      static constexpr bool b1 = ((requires(const EvaluatorType& e) {
                                    e(std::declval<element_index<Space>>());
                                  }) &&
                                  (ElementSpaceConcept<Space>));
      static constexpr bool b2 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<element_workspace<Space>>(),
               std::declval<element_index<Space>>());
           }) &&
           (ElementSpaceConcept<Space> && hasElementWorkspace<Space>));
      static constexpr bool b3 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<cell_index<Space>>(),
               std::declval<quadrature_point_index<Space>>());
           }) &&
           (QuadratureSpaceConcept<Space>));
      static constexpr bool b4 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<cell_workspace<Space>>(),
               std::declval<cell_index<Space>>(),
               std::declval<quadrature_point_index<Space>>());
           }) &&
           (QuadratureSpaceConcept<Space> && hasCellWorkspace<Space>));
      using ResultType1 =
          typename EvaluatorResultQueryImplementation1<b1, EvaluatorType>::type;
      using ResultType2 =
          typename EvaluatorResultQueryImplementation2<b2, EvaluatorType>::type;
      using ResultType3 =
          typename EvaluatorResultQueryImplementation3<b3, EvaluatorType>::type;
      using ResultType4 =
          typename EvaluatorResultQueryImplementation4<b4, EvaluatorType>::type;
      using type = std::conditional_t<
          b1,
          ResultType1,
          std::conditional_t<b2,
                             ResultType2,
                             std::conditional_t<b3, ResultType3, ResultType4>>>;
    };

    template <typename EvaluatorType>
        struct EvaluatorResultQuery
        : EvaluatorResultQueryImplementation < requires(EvaluatorType& e) {
      getSpace(e);
    }, EvaluatorType > {};

  }  // namespace internals

  /*!
   * \brief a concept that must satisfy an evaluator
   */
  template <typename EvaluatorType>
  concept EvaluatorConcept = std::is_move_constructible_v<EvaluatorType> &&
      std::is_copy_constructible_v<EvaluatorType> && SpaceConcept<
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>> &&
      requires(const EvaluatorType& e) {
    { getNumberOfComponents(e) } -> std::same_as<mgis::size_type>;
  } && requires(const EvaluatorType& e, AbstractErrorHandler& eh) {
    { check(eh, e) } -> std::same_as<bool>;
  } &&((internals::EvaluatorResultQuery<EvaluatorType>::b1) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b2) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b3) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b4)) &&
      (internals::EvaluatorResultQuery<EvaluatorType>::b1
           ? std::same_as<
                 typename internals::EvaluatorResultQuery<
                     EvaluatorType>::ResultType1,
                 typename internals::EvaluatorResultQuery<EvaluatorType>::type>
           : true) &&
      (internals::EvaluatorResultQuery<EvaluatorType>::b2
           ? std::same_as<
                 typename internals::EvaluatorResultQuery<
                     EvaluatorType>::ResultType2,
                 typename internals::EvaluatorResultQuery<EvaluatorType>::type>
           : true) &&
      (internals::EvaluatorResultQuery<EvaluatorType>::b3
           ? std::same_as<
                 typename internals::EvaluatorResultQuery<
                     EvaluatorType>::ResultType3,
                 typename internals::EvaluatorResultQuery<EvaluatorType>::type>
           : true) &&
      (internals::EvaluatorResultQuery<EvaluatorType>::b4
           ? std::same_as<
                 typename internals::EvaluatorResultQuery<
                     EvaluatorType>::ResultType4,
                 typename internals::EvaluatorResultQuery<EvaluatorType>::type>
           : true);

  //! concept defining evaluators working on an element space
  template <typename EvaluatorType>
  concept ElementEvaluatorConcept = EvaluatorConcept<EvaluatorType> &&
      ((internals::EvaluatorResultQuery<EvaluatorType>::b1) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b2));

  //! concept defining evaluators working on a quadrature space
  template <typename EvaluatorType>
  concept QuadratureEvaluatorConcept = EvaluatorConcept<EvaluatorType> &&
      ((internals::EvaluatorResultQuery<EvaluatorType>::b3) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b4));

  //
  template <EvaluatorConcept EvaluatorType>
  using evaluator_space =
      std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;

  /*!
   * \brief type of the result of an evaluator
   */
  template <EvaluatorConcept EvaluatorType>
  using evaluator_result =
      typename internals::EvaluatorResultQuery<EvaluatorType>::type;

  namespace internals {

    template <typename T>
    concept is_pointer_to_real = (std::same_as<real*, T>) ||
                                 (std::same_as<const real*, T>);

  }  // end of   namespace internals

  template <EvaluatorConcept EvaluatorType>
  inline constexpr auto isEvaluatorResultTypeMappable =
      requires(evaluator_result<EvaluatorType> rf) {
    { rf.data() } -> internals::is_pointer_to_real;
  };

  namespace internals {

    /*!
     * This helper function allows to disambiguate the call to
     * the getSpace function
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr decltype(auto) disambiguateGetSpace(const EvaluatorType&);
    /*!
     * This helper function allows to disambiguate the call to
     * the check function
     */
    template <EvaluatorConcept EvaluatorType>
    [[nodiscard]] constexpr bool disambiguateCheck(AbstractErrorHandler&,
                                                   const EvaluatorType&);
    /*!
     * This helper function allows to disambiguate the call to
     * the getNumberOfComponents function
     */
    template <EvaluatorConcept EvaluatorType>
    [[nodiscard]] constexpr mgis::size_type disambiguateGetNumberOfComponents(
        EvaluatorType&);

    /*!
     * \brief number of components of a type when known at compile-time,
     * dynamic_extent otherwise
     */
    template <EvaluatorConcept EvaluatorType>
    struct NumberOfComponents<EvaluatorType> {
      static constexpr size_type value = internals::CompileTimeSize<
          std::decay_t<evaluator_result<EvaluatorType>>>::value;
    };

  }  // namespace internals

}  // end of namespace mgis::function

#include "MGIS/Function/EvaluatorConcept.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX */
