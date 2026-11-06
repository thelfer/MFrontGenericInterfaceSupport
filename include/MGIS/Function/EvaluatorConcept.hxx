/*!
 * \file   MGIS/Function/EvaluatorConcept.hxx
 * \brief File which defines the concept an evaluator must satisfy and the
 * helpers used to query the space on which an evaluator is defined and the type
 * of the values it returns.
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

#ifndef LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX
#define LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX

#include <concepts>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function {

  namespace internals {

    /*!
     * \brief helper class used to determine the type returned by an evaluator
     * called with an element index, which is an invalid type if the evaluator
     * cannot be called this way
     *
     * \tparam b: false if the evaluator cannot be called with an element index
     * \tparam EvaluatorType: evaluator type
     */
    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation1 {
      //! \brief type returned when the evaluator cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the evaluator
      using type = InvalidResult;
    };

    /*!
     * \brief specialization used when the evaluator can be called this way
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation1<true, EvaluatorType> {
      //! \brief space on which the evaluator is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      //! \brief type of the value returned by the evaluator
      using type =
          std::invoke_result_t<const EvaluatorType, element_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by an evaluator
     * called with an element workspace and an element index, which is an
     * invalid type if the evaluator cannot be called this way
     *
     * \tparam b: false if the evaluator cannot be called this way
     * \tparam EvaluatorType: evaluator type
     */
    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation2 {
      //! \brief type returned when the evaluator cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the evaluator
      using type = InvalidResult;
    };

    /*!
     * \brief specialization used when the evaluator can be called this way
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation2<true, EvaluatorType> {
      //! \brief space on which the evaluator is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      //! \brief type of the value returned by the evaluator
      using type = std::invoke_result_t<const EvaluatorType,
                                        element_workspace<Space>,
                                        element_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by an evaluator
     * called with a cell index and a quadrature point index, which is an
     * invalid type if the evaluator cannot be called this way
     *
     * \tparam b: false if the evaluator cannot be called this way
     * \tparam EvaluatorType: evaluator type
     */
    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation3 {
      //! \brief type returned when the evaluator cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the evaluator
      using type = InvalidResult;
    };

    /*!
     * \brief specialization used when the evaluator can be called this way
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation3<true, EvaluatorType> {
      //! \brief space on which the evaluator is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      //! \brief type of the value returned by the evaluator
      using type = std::invoke_result_t<const EvaluatorType,
                                        cell_index<Space>,
                                        quadrature_point_index<Space>>;
    };

    /*!
     * \brief helper class used to determine the type returned by an evaluator
     * called with a cell workspace, a cell index and a quadrature point index,
     * which is an invalid type if the evaluator cannot be called this way
     *
     * \tparam b: false if the evaluator cannot be called this way
     * \tparam EvaluatorType: evaluator type
     */
    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation4 {
      //! \brief type returned when the evaluator cannot be called this way
      struct InvalidResult;
      //! \brief type of the value returned by the evaluator
      using type = InvalidResult;
    };

    /*!
     * \brief specialization used when the evaluator can be called this way
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation4<true, EvaluatorType> {
      //! \brief space on which the evaluator is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      //! \brief type of the value returned by the evaluator
      using type = std::invoke_result_t<const EvaluatorType,
                                        cell_workspace<Space>,
                                        cell_index<Space>,
                                        quadrature_point_index<Space>>;
    };

    /*!
     * \brief helper class used to determine how an evaluator can be called and
     * the type of the values it returns
     *
     * \tparam b: false if the evaluator does not provide a `getSpace` function
     * \tparam EvaluatorType: evaluator type
     */
    template <bool, typename EvaluatorType>
    struct EvaluatorResultQueryImplementation {
      //! \brief type used for the queries which failed
      struct InvalidResult;
      //! \brief the evaluator can be called with an element index
      static constexpr bool b1 = false;
      //! \brief the evaluator can be called with an element workspace and an
      /// element index
      static constexpr bool b2 = false;
      //! \brief the evaluator can be called with a cell index and a quadrature
      /// point index
      static constexpr bool b3 = false;
      //! \brief the evaluator can be called with a cell workspace, a cell index
      /// and a quadrature point index
      static constexpr bool b4 = false;
      //! \brief type returned by an evaluator called with an element index
      using ResultType1 = InvalidResult;
      //! \brief type returned by an evaluator called with an element workspace
      /// and an element index
      using ResultType2 = InvalidResult;
      //! \brief type returned by an evaluator called with a cell index and a
      /// quadrature point index
      using ResultType3 = InvalidResult;
      //! \brief type returned by an evaluator called with a cell workspace, a
      /// cell index and a quadrature point index
      using ResultType4 = InvalidResult;
    };

    /*!
     * \brief specialization used when the evaluator provides a `getSpace`
     * function
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
    struct EvaluatorResultQueryImplementation<true, EvaluatorType> {
      //! \brief space on which the evaluator is defined
      using Space =
          std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;
      //! \brief the evaluator can be called with an element index
      static constexpr bool b1 = ((requires(const EvaluatorType& e) {
                                    e(std::declval<element_index<Space>>());
                                  }) &&
                                  (ElementSpaceConcept<Space>));
      //! \brief the evaluator can be called with an element workspace and an
      /// element index
      static constexpr bool b2 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<element_workspace<Space>>(),
               std::declval<element_index<Space>>());
           }) &&
           (ElementSpaceConcept<Space> && hasElementWorkspace<Space>));
      //! \brief the evaluator can be called with a cell index and a quadrature
      /// point index
      static constexpr bool b3 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<cell_index<Space>>(),
               std::declval<quadrature_point_index<Space>>());
           }) &&
           (QuadratureSpaceConcept<Space>));
      //! \brief the evaluator can be called with a cell workspace, a cell index
      /// and a quadrature point index
      static constexpr bool b4 =
          ((requires(const EvaluatorType& e) {
             e(std::declval<cell_workspace<Space>>(),
               std::declval<cell_index<Space>>(),
               std::declval<quadrature_point_index<Space>>());
           }) &&
           (QuadratureSpaceConcept<Space> && hasCellWorkspace<Space>));
      //! \brief type returned by an evaluator called with an element index
      using ResultType1 =
          typename EvaluatorResultQueryImplementation1<b1, EvaluatorType>::type;
      //! \brief type returned by an evaluator called with an element workspace
      /// and an element index
      using ResultType2 =
          typename EvaluatorResultQueryImplementation2<b2, EvaluatorType>::type;
      //! \brief type returned by an evaluator called with a cell index and a
      /// quadrature point index
      using ResultType3 =
          typename EvaluatorResultQueryImplementation3<b3, EvaluatorType>::type;
      //! \brief type returned by an evaluator called with a cell workspace, a
      /// cell index and a quadrature point index
      using ResultType4 =
          typename EvaluatorResultQueryImplementation4<b4, EvaluatorType>::type;
      //! \brief type of the values returned by the evaluator, which is the type
      /// returned by the first valid call
      using type = std::conditional_t<
          b1,
          ResultType1,
          std::conditional_t<b2,
                             ResultType2,
                             std::conditional_t<b3, ResultType3, ResultType4>>>;
    };

    /*!
     * \brief alias used to determine how an evaluator can be called and the
     * type of the values it returns
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <typename EvaluatorType>
        struct EvaluatorResultQuery
        : EvaluatorResultQueryImplementation < requires(EvaluatorType& e) {
      getSpace(e);
    }, EvaluatorType > {};

  }  // namespace internals

  /*!
   * \brief a concept that must satisfy an evaluator
   *
   * An evaluator is an object which:
   * - can be moved and copied,
   * - is defined on a space, i.e. for which the `getSpace` function is
   *   available,
   * - provides the number of components of the values it returns through the
   *   `getNumberOfComponents` function,
   * - can check the values it computes through the `check` function,
   * - can be called with either an element index, or an element workspace and
   *   an element index, or a cell index and a quadrature point index, or a cell
   *   workspace, a cell index and a quadrature point index. All the valid calls
   *   return the same type.
   *
   * Contrary to a function, an evaluator is called with a constant object and
   * returns constant values.
   *
   * \tparam EvaluatorType: evaluator type
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

  /*!
   * \brief concept defining evaluators working on an element space
   *
   * \tparam EvaluatorType: evaluator type
   */
  template <typename EvaluatorType>
  concept ElementEvaluatorConcept = EvaluatorConcept<EvaluatorType> &&
      ((internals::EvaluatorResultQuery<EvaluatorType>::b1) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b2));

  /*!
   * \brief concept defining evaluators working on a quadrature space
   *
   * \tparam EvaluatorType: evaluator type
   */
  template <typename EvaluatorType>
  concept QuadratureEvaluatorConcept = EvaluatorConcept<EvaluatorType> &&
      ((internals::EvaluatorResultQuery<EvaluatorType>::b3) ||
       (internals::EvaluatorResultQuery<EvaluatorType>::b4));

  /*!
   * \brief alias giving the type of the space on which an evaluator is defined
   *
   * \tparam EvaluatorType: evaluator type
   */
  template <EvaluatorConcept EvaluatorType>
  using evaluator_space =
      std::decay_t<decltype(getSpace(std::declval<EvaluatorType>()))>;

  /*!
   * \brief type of the result of an evaluator, i.e. the type of the constant
   * values returned by the evaluator, which can be a scalar value, a constant
   * `std::span` or a constant `std::array`
   *
   * \tparam EvaluatorType: evaluator type
   */
  template <EvaluatorConcept EvaluatorType>
  using evaluator_result =
      typename internals::EvaluatorResultQuery<EvaluatorType>::type;

  /*!
   * \brief check if the given evaluators share the same space
   *
   * \param[in] ctx: error handler used to report the errors
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   * \return true if the two evaluators share the same space
   */
  [[nodiscard]] constexpr bool checkMatchingSpaces(
      AbstractErrorHandler& ctx,
      const EvaluatorConcept auto& e1,
      const EvaluatorConcept auto& e2);

  namespace internals {

    /*!
     * \brief concept satisfied by a pointer to a scalar value, possibly
     * constant
     *
     * \tparam T: type
     */
    template <typename T>
    concept is_pointer_to_real = (std::same_as<real*, T>) ||
                                 (std::same_as<const real*, T>);

    /*!
     * \brief concept satisfied by `real`, possibly const or a reference
     *
     * \tparam T: type to be tested
     */
    template <typename T>
    concept is_real = std::same_as<std::remove_cvref_t<T>, real>;

  }  // end of   namespace internals

  /*!
   * \brief returns a pointer to the values held by an object with a `data`
   * method, such as `std::span`, `std::array` or the tensorial objects of
   * `TFEL/Math`
   *
   * \tparam ResultType: type of the object
   * \param[in] r: object exposing a `data` method returning a pointer to
   * scalar values
   * \return a pointer to the values held by the object
   *
   * \note `getDataPointer` is used by the views and modifiers to map the
   * results of functions and evaluators. It can be overloaded to support other
   * result types.
   */
  template <typename ResultType>
  [[nodiscard]] constexpr auto getDataPointer(ResultType&& r) requires(
      requires(ResultType & rr) {
        { rr.data() } -> internals::is_pointer_to_real;
      }) {
    return r.data();
  }
  /*!
   * \brief returns a pointer to a scalar value
   *
   * \tparam ScalarType: type of the scalar value, which may be a constant
   * reference or a mutable reference
   * \param[in] r: scalar value
   * \return a pointer to the scalar value
   */
  template <internals::is_real ScalarType>
  [[nodiscard]] constexpr auto getDataPointer(ScalarType&& r) MGIS_NOEXCEPT {
    return &r;
  }

  /*!
   * \brief true if the result of an evaluator is a mappable type, i.e. a type
   * for which the `getDataPointer` function returns a pointer to the values it
   * holds
   *
   * \tparam EvaluatorType: evaluator type
   */
  template <EvaluatorConcept EvaluatorType>
  inline constexpr auto isEvaluatorResultTypeMappable =
      requires(evaluator_result<EvaluatorType> rf) {
    { getDataPointer(rf) } -> internals::is_pointer_to_real;
  };

  namespace internals {

    /*!
     * \brief This helper function allows to disambiguate the call to
     * the getSpace function
     *
     * \tparam EvaluatorType: evaluator type
     * \param[in] e: evaluator
     * \return the space on which the evaluator is defined
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr decltype(auto) disambiguateGetSpace(const EvaluatorType& e);
    /*!
     * \brief This helper function allows to disambiguate the call to
     * the check function
     *
     * \tparam EvaluatorType: evaluator type
     * \param[in] eh: error handler
     * \param[in] e: evaluator
     * \return true if the values computed by the evaluator are valid
     */
    template <EvaluatorConcept EvaluatorType>
    [[nodiscard]] constexpr bool disambiguateCheck(AbstractErrorHandler& eh,
                                                   const EvaluatorType& e);
    /*!
     * \brief This helper function allows to disambiguate the call to
     * the getNumberOfComponents function
     *
     * \tparam EvaluatorType: evaluator type
     * \param[in] e: evaluator
     * \return the number of components of the values returned by the evaluator
     */
    template <EvaluatorConcept EvaluatorType>
    [[nodiscard]] constexpr mgis::size_type disambiguateGetNumberOfComponents(
        const EvaluatorType& e);

    /*!
     * \brief number of components of an evaluator when known at compile-time,
     * dynamic_extent otherwise
     *
     * \tparam EvaluatorType: evaluator type
     */
    template <EvaluatorConcept EvaluatorType>
    struct NumberOfComponents<EvaluatorType> {
      //! \brief number of components of the values returned by the evaluator
      static constexpr size_type value = internals::CompileTimeSize<
          std::decay_t<evaluator_result<EvaluatorType>>>::value;
    };

    /*!
     * \brief checks if the number of components of an evaluator is compatible
     * with the number of components expected by a modifier
     *
     * The check is always satisfied if the number of components of the
     * evaluator is only known at run time.
     *
     * \tparam EvaluatorType: evaluator type
     * \tparam N: number of components expected by the modifier
     * \return true if the number of components is compatible
     */
    template <EvaluatorConcept EvaluatorType, size_type N>
    [[nodiscard]] constexpr bool checkEvaluatorNumberOfComponentsCompatibility()
        MGIS_NOEXCEPT {
      constexpr auto Nc = number_of_components<EvaluatorType>;
      if constexpr (Nc == dynamic_extent) {
        return true;
      } else {
        return N == Nc;
      }
    }  // end of checkEvaluatorNumberOfComponentsCompatibility

  }  // namespace internals

}  // end of namespace mgis::function

#include "MGIS/Function/EvaluatorConcept.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORCONCEPT_HXX */
