/*!
 * \file   MGIS/Function/TFEL/QuantityModifier.hxx
 * \brief  File which defines the `QuantityModifier` class, an evaluator
 * returning the values computed by an evaluator with one component as TFEL's
 * quantities, i.e. scalar values associated with a physical unit.
 * \author Thomas Helfer
 * \date   16/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_HXX

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/TFEL/General.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values computed by an evaluator as
   * quantities
   *
   * A quantity is a scalar value associated with a physical unit. The modifier
   * is built from an evaluator with exactly one component. Its call operators
   * return the value computed by the evaluator at a given integration point:
   *
   * - a read-only quantity (a `const_qt_ref` object of `TFEL/Math`)
   *   referencing the value computed by the evaluator when this value is
   *   stored in memory owned by another object,
   * - a quantity owning its value (a `qt` object of `TFEL/Math`), in which
   *   case the value computed by the evaluator is copied.
   *
   * Only the call operators which match the way the underlying evaluator can
   * be called are available.
   *
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam UnitType: unit of the quantities returned by the modifier
   */
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct QuantityModifier : private PreconditionsChecker<
                                QuantityModifier<EvaluatorType, UnitType>> {
    //! \brief type of the space on which the underlying evaluator is defined
    using Space = evaluator_space<EvaluatorType>;
    /*!
     * \brief method checking that the preconditions of the constructor are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] values: evaluator used to build the modifier
     * \return true if the number of components of the evaluator is equal to one
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const EvaluatorType&);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] values: evaluator used to build the modifier
     */
    constexpr QuantityModifier(const EvaluatorType&);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: evaluator used to build the modifier
     */
    template <bool doPreconditionsCheck>
    constexpr QuantityModifier(const PreconditionsCheck<doPreconditionsCheck>&,
                               const EvaluatorType&);
    /*!
     * \brief perform consistency checks
     * \param[in] ctx: error handler used to report the errors
     * \return true if the values computed by the underlying evaluator are
     * valid and if its number of components is equal to one
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying evaluator is defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \return the number of components of the values returned by the
     * modifier, which is one
     */
    [[nodiscard]] constexpr size_type getNumberOfComponents() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     * \return the value computed by the evaluator at the given integration
     * point as a quantity
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
     * \return the value computed by the evaluator at the given integration
     * point as a quantity
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
     * \return the value computed by the evaluator at the given integration
     * point as a quantity
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
     * \return the value computed by the evaluator at the given integration
     * point as a quantity
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
  };  // end of QuantityModifier

  /*!
   * \brief return the underlying space
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam UnitType: unit of the quantities returned by the modifier
   * \param[in] e: quantity modifier
   * \return the space on which the modifier is defined
   */
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const QuantityModifier<EvaluatorType, UnitType>&);
  /*!
   * \brief perform consistency checks
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam UnitType: unit of the quantities returned by the modifier
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: quantity modifier
   * \return true if the values computed by the modifier are valid
   */
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const QuantityModifier<EvaluatorType, UnitType>&);
  /*!
   * \brief return the number of components
   * \tparam EvaluatorType: type of the underlying evaluator
   * \tparam UnitType: unit of the quantities returned by the modifier
   * \param[in] e: quantity modifier
   * \return the number of components of the modifier, which is one
   */
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const QuantityModifier<EvaluatorType, UnitType>&);

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/QuantityModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_HXX */
