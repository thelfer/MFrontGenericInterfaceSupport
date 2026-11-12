/*!
 * \file   MGIS/Function/FixedSizeModifier.hxx
 * \brief File which defines the `FixedSizeModifier` class, a modifier returning
 * the values computed by an evaluator with a fixed number of components as a
 * fixed size span or as a scalar.
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

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX

#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of a function view as a fixed
   * size span or a scalar
   *
   * The modifier is built from an evaluator with exactly `N` components. Its
   * call operators return the values computed by the evaluator at a given
   * integration point:
   *
   * - a scalar if `N` is equal to 1,
   * - a `std::span<const real, N>` if `N` is greater than 1 and the values
   *   returned by the evaluator are stored in contiguous memory,
   * - a `std::array<real, N>` otherwise, in which case the values are copied.
   *
   * Only the call operators which match the way the underlying evaluator can
   * be called are available.
   *
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam N: number of components of the values returned by the modifier,
   * which must be strictly positive
   */
  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) struct FixedSizeModifier
      : private PreconditionsChecker<FixedSizeModifier<EvaluatorType, N>> {
    //! \brief type of the space on which the underlying evaluator is defined
    using Space = evaluator_space<EvaluatorType>;
    /*!
     * \brief method checking that the preconditions of the constructor are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] e: evaluator used to build the modifier
     * \return true if the number of components of the evaluator is equal to
     * `N`
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const EvaluatorType&);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] values: evaluator used to build the modifier
     */
    constexpr FixedSizeModifier(const EvaluatorType&);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: evaluator used to build the modifier
     */
    template <bool doPreconditionsCheck>
    constexpr FixedSizeModifier(const PreconditionsCheck<doPreconditionsCheck>&,
                                const EvaluatorType&);
    /*!
     * \brief perform consistency checks
     * \param[in] ctx: error handler used to report the errors
     * \return true if the values computed by the underlying evaluator are
     * valid and if its number of components is equal to `N`
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying evaluator is defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \return the number of components of the values returned by the
     * modifier, which is `N`
     */
    [[nodiscard]] constexpr size_type getNumberOfComponents() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     * \return the values computed by the evaluator at the given integration
     * point: a scalar if `N` is equal to 1, fixed size read-only values
     * otherwise
     * \note only available if the underlying evaluator can be called with an
     * element index and if the values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b1) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the values computed by the evaluator at the given integration
     * point: a scalar if `N` is equal to 1, fixed size read-only values
     * otherwise
     * \note only available if the underlying evaluator can be called with an
     * element workspace and an element index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b2) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values computed by the evaluator at the given integration
     * point: a scalar if `N` is equal to 1, fixed size read-only values
     * otherwise
     * \note only available if the underlying evaluator can be called with a
     * cell index and a quadrature point index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b3) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values computed by the evaluator at the given integration
     * point: a scalar if `N` is equal to 1, fixed size read-only values
     * otherwise
     * \note only available if the underlying evaluator can be called with a
     * cell workspace, a cell index and a quadrature point index and if the
     * values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b4) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));

   private:
    //! \brief underlying evaluator
    EvaluatorType evaluator;
  };  // end of FixedSizeModifier

  /*!
   * \brief create a fixed size view from an evaluator
   * \tparam N: number of components of the view
   * \tparam EvaluatorType: type of the evaluator, which may be a qualified type
   * \param[in] f: evaluator
   * \return a fixed size modifier of the evaluator
   */
  template <size_type N, typename EvaluatorType>
  [[nodiscard]] constexpr auto view(const EvaluatorType&) requires(
      (N > 0) && (N != dynamic_extent) &&  //
      (EvaluatorConcept<std::decay_t<EvaluatorType>>));

  /*!
   * \brief return the underlying space
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam N: number of components of the modifier
   * \param[in] e: fixed size modifier
   * \return the space on which the modifier is defined
   */
  template <EvaluatorConcept EvaluatorType, size_type N>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const FixedSizeModifier<EvaluatorType, N>&);
  /*!
   * \brief perform consistency checks
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam N: number of components of the modifier
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: fixed size modifier
   * \return true if the values computed by the modifier are valid
   */
  template <EvaluatorConcept EvaluatorType, size_type N>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const FixedSizeModifier<EvaluatorType, N>&);
  /*!
   * \brief return the number of components
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam N: number of components of the modifier
   * \param[in] e: fixed size modifier
   * \return the number of components of the modifier, which is `N`
   */
  template <EvaluatorConcept EvaluatorType, size_type N>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const FixedSizeModifier<EvaluatorType, N>&);

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX */
