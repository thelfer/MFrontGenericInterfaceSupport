/*!
 * \file   MGIS/Function/FixedSizeView.hxx
 * \brief File which defines the `FixedSizeView` class, a view returning the
 * values of a function with a fixed number of components as a fixed size span
 * or as a scalar.
 * \author Thomas Helfer
 * \date   07/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX

#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*!
   * \brief a modifier returning the values of a function view as a fixed size
   * span or a scalar
   *
   * The view is built from a function with exactly `N` components. Its call
   * operators return the values associated with an integration point:
   *
   * - a scalar if `N` is equal to 1 for the constant call operators and a
   *   mutable reference to a scalar for the non-constant ones,
   * - a `std::span<const real, N>` for the constant call operators and a
   *   `std::span<real, N>` for the non-constant ones if `N` is greater
   *   than 1.
   *
   * Only the call operators which match the way the underlying function can be
   * called are available.
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam N: number of components of the values returned by the view, which
   * must be strictly positive
   */
  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0) struct FixedSizeView
      : private PreconditionsChecker<FixedSizeView<FunctionType, N>> {
    //! \brief type of the space on which the underlying function is defined
    using Space = function_space<FunctionType>;
    /*!
     * \brief value returned by the non-const call operators: a mutable
     * reference to a scalar if `N` is equal to 1, a mutable fixed size span
     * otherwise
     */
    using mutable_value_type =
        std::conditional_t<N == 1, real&, std::span<real, N>>;
    /*!
     * \brief a simple alias used to workaround what seems to be a bug in
     * gcc 16.x: the underlying function is passed to the constructors by value
     * if it is a lightweight view and by reference otherwise
     */
    using ConstructorArgumentType =
        std::conditional_t<LightweightViewConcept<FunctionType>,
                           FunctionType,
                           FunctionType&>;
    /*!
     * \brief method checking that the preconditions of the constructor are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] values: function used to build the view
     * \return true if the number of components of the function is equal to `N`
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const FunctionType&);

    /*!
     * \brief constructor, checking the preconditions
     * \param[in] values: function used to build the view
     */
    constexpr FixedSizeView(ConstructorArgumentType);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function used to build the view
     */
    template <bool doPreconditionsCheck>
    constexpr FixedSizeView(const PreconditionsCheck<doPreconditionsCheck>&,
                            ConstructorArgumentType);
    /*!
     * \brief perform consistency checks
     * \param[in] eh: error handler used to report the errors
     * \return true if the number of components of the underlying function is
     * equal to `N`
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying function is defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \return the number of components of the values returned by the view,
     * which is `N`
     */
    [[nodiscard]] constexpr size_type getNumberOfComponents() const
        MGIS_NOEXCEPT;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point: a
     * scalar if `N` is equal to 1, a `std::span<const real, N>` otherwise
     * \note only available if the underlying function can be called with an
     * element index and if the values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point: a
     * scalar if `N` is equal to 1, a `std::span<const real, N>` otherwise
     * \note only available if the underlying function can be called with an
     * element workspace and an element index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the value of the function at the given integration point: a
     * scalar if `N` is equal to 1, a `std::span<const real, N>` otherwise
     * \note only available if the underlying function can be called with a
     * cell index and a quadrature point index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the value of the function at the given integration point: a
     * scalar if `N` is equal to 1, a `std::span<const real, N>` otherwise
     * \note only available if the underlying function can be called with a
     * cell workspace, a cell index and a quadrature point index and if the
     * values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point, which
     * may be modified: a mutable reference to a scalar if `N` is equal to 1, a
     * `std::span<real, N>` otherwise
     * \note only available if the underlying function can be called with an
     * element index and if the values it returns can be mapped to mutable
     * memory
     */
    [[nodiscard]] constexpr mutable_value_type operator()(
        const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point, which
     * may be modified: a mutable reference to a scalar if `N` is equal to 1, a
     * `std::span<real, N>` otherwise
     * \note only available if the underlying function can be called with an
     * element workspace and an element index and if the values it returns can
     * be mapped to mutable memory
     */
    [[nodiscard]] constexpr mutable_value_type operator()(
        const element_workspace<Space>&,
        const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the value of the function at the given integration point, which
     * may be modified: a mutable reference to a scalar if `N` is equal to 1, a
     * `std::span<real, N>` otherwise
     * \note only available if the underlying function can be called with a
     * cell index and a quadrature point index and if the values it returns can
     * be mapped to mutable memory
     */
    [[nodiscard]] constexpr mutable_value_type operator()(
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the value of the function at the given integration point, which
     * may be modified: a mutable reference to a scalar if `N` is equal to 1, a
     * `std::span<real, N>` otherwise
     * \note only available if the underlying function can be called with a
     * cell workspace, a cell index and a quadrature point index and if the
     * values it returns can be mapped to mutable memory
     */
    [[nodiscard]] constexpr mutable_value_type operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionResultTypeMappable<FunctionType>));

   private:
    //! \brief view of the underlying function
    function_view<FunctionType> function;
  };  // end of FixedSizeView

  /*!
   * \brief specialisation stating that a `FixedSizeView` is a lightweight
   * view, i.e. that a copy of the view is cheap and can be used to evaluate
   * the function
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam N: number of components of the view
   */
  template <FunctionConcept FunctionType, size_type N>
  struct LightweightViewTraits<FixedSizeView<FunctionType, N>>
      : std::true_type {};

  /*!
   * \brief convert a function to a fixed size view
   * \tparam N: number of components of the view
   * \tparam FunctionType: type of the function, which may be a qualified type
   * \param[in] f: function
   * \return a fixed size view of the function
   * \note this overload applies to non-constant functions and to functions
   * which are their own view
   */
  template <size_type N, typename FunctionType>
  [[nodiscard]] constexpr auto view(FunctionType&) requires(
      (N > 0) && (N != dynamic_extent) &&              //
      (FunctionConcept<std::decay_t<FunctionType>>)&&  //
      ((!std::is_const_v<FunctionType>) ||
       (LightweightViewConcept<std::decay_t<FunctionType>>)));
  /*!
   * \brief convert a const function to an immutable fixed size view
   * \tparam N: number of components of the view
   * \tparam FunctionType: const type of the function
   * \param[in] f: function
   * \return a fixed size view of a view of the function, which only returns
   * constant values
   * \note this overload applies to constant functions which are neither
   * evaluators nor lightweight views and for which a view can be created
   */
  template <size_type N, typename FunctionType>
  [[nodiscard]] constexpr auto view(FunctionType&) requires(
      (N > 0) && (N != dynamic_extent) &&                //
      (FunctionConcept<std::decay_t<FunctionType>>)&&    //
      (!EvaluatorConcept<std::decay_t<FunctionType>>)&&  //
      (!LightweightViewConcept<std::decay_t<FunctionType>>)&&  //
      (std::is_const_v<FunctionType>)&&                  //
      (requires(FunctionType & rf) { view(rf); }));

  /*!
   * \brief return the underlying space
   * \tparam FunctionType: type of the underlying function
   * \tparam N: number of components of the view
   * \param[in] v: fixed size view
   * \return the space on which the view is defined
   */
  template <FunctionConcept FunctionType, size_type N>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const FixedSizeView<FunctionType, N>&);
  /*!
   * \brief perform consistency checks
   * \tparam FunctionType: type of the underlying function
   * \tparam N: number of components of the view
   * \param[in] eh: error handler used to report the errors
   * \param[in] v: fixed size view
   * \return true if the number of components of the view is equal to `N`
   */
  template <FunctionConcept FunctionType, size_type N>
  [[nodiscard]] constexpr bool check(AbstractErrorHandler&,
                                     const FixedSizeView<FunctionType, N>&);
  /*!
   * \brief return the number of components
   * \tparam FunctionType: type of the underlying function
   * \tparam N: number of components of the view
   * \param[in] v: fixed size view
   * \return the number of components of the view, which is `N`
   */
  template <FunctionConcept FunctionType, size_type N>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const FixedSizeView<FunctionType, N>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeView.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX */
