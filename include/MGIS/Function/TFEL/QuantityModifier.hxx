/*!
 * \file   MGIS/Function/QuantityModifier.hxx
 * \brief
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

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of a
   * function view as a fixed size span or a scalar
   *
   * \tparam Space: functional space
   * \tparam N: size of the returned value
   */
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct QuantityModifier : private PreconditionsChecker<
                                QuantityModifier<EvaluatorType, UnitType>> {
    //
    using Space = evaluator_space<EvaluatorType>;
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] eh: error handler
     * \param[in] e: evaluator
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const EvaluatorType&);
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    constexpr QuantityModifier(const EvaluatorType&);
    /*!
     * \brief constructor
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function
     */
    template <bool doPreconditionsCheck>
    constexpr QuantityModifier(const PreconditionsCheck<doPreconditionsCheck>&,
                               const EvaluatorType&);
    //! \brief perform consistency checks
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    //! \brief return the underlying  space
    [[nodiscard]] decltype(auto) getSpace() const;
    //! \return the number of components
    [[nodiscard]] constexpr size_type getNumberOfComponents() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b1) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b2) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b3) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires((internals::EvaluatorResultQuery<EvaluatorType>::b4) &&
                 (isEvaluatorResultTypeMappable<EvaluatorType>));

   private:
    //! \brief underlying function
    EvaluatorType evaluator;
  };  // end of QuantityModifier

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] decltype(auto) getSpace(
      const QuantityModifier<EvaluatorType, UnitType>&);
  //! \brief perform consistency checks
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const QuantityModifier<EvaluatorType, UnitType>&);
  //! \return the number of components
  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const QuantityModifier<EvaluatorType, UnitType>&);

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/QuantityModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_HXX */
