/*!
 * \file   MGIS/Function/TFEL/TensorModifier.hxx
 * \brief  File which defines the `TensorModifier` class, an evaluator
 * returning the values computed by an evaluator as TFEL's tensorial objects.
 * \author Thomas Helfer
 * \date   13/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/EvaluatorModifierBase.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a variable stating if the requirements of the `TensorModifier`
   * class are met
   *
   * The values returned by the evaluator must be mappable to memory and their
   * number of components must be compatible with the number of components of
   * the tensorial type, i.e. be equal to it or be only known at run time.
   *
   * \tparam TensorType: type of the tensorial object to be exposed
   * \tparam EvaluatorType: type of the underlying evaluator
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  inline constexpr bool areTensorModifierRequirementsSatisfied =
      ((isEvaluatorResultTypeMappable<EvaluatorType>)&&  //
       (internals::checkEvaluatorNumberOfComponentsCompatibility<
           EvaluatorType,
           compile_time_size<TensorType>>()));

  /*!
   * \brief an evaluator returning the values computed by an evaluator as TFEL's
   * tensorial objects
   *
   * The modifier is built from an evaluator whose number of components must be
   * compatible with the tensorial type `TensorType`. Its call operators pass
   * the value computed by the evaluator at a given integration point to the
   * `apply` member function, which returns it as a tensorial object: a
   * read-only tensor viewing the computed values if they are stored in memory
   * owned by another object, a tensor owning a copy of these values otherwise.
   *
   * Only the call operators which match the way the underlying evaluator can be
   * called are available.
   *
   * \tparam TensorType: type of the tensorial object to be exposed, which must
   * match the `TensorConcept` concept
   * \tparam EvaluatorType: type of the underlying evaluator
   *
   * \note only available if the requirements stated by the
   * `areTensorModifierRequirementsSatisfied` variable are satisfied
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      struct TensorModifier
      : EvaluatorModifierBase<TensorModifier<TensorType, EvaluatorType>,
                              EvaluatorType> {
    /*!
     * \brief constructors inherited from the `EvaluatorModifierBase` base class
     */
    using EvaluatorModifierBase<TensorModifier,
                                EvaluatorType>::EvaluatorModifierBase;
    /*!
     * \brief perform consistency checks
     * \param[in] eh: error handler the detected errors are reported to
     * \return true if the values computed by the underlying evaluator are valid
     * and if its number of components is compatible with the tensorial type
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \return the number of components of the values returned by the modifier,
     * which is the number of components of the tensorial type
     */
    constexpr size_type getNumberOfComponents() const MGIS_NOEXCEPT;
    /*!
     * \brief transform the value computed by the underlying evaluator at a
     * given integration point into a tensorial object
     * \param[in] values: value computed by the underlying evaluator at the
     * given integration point
     * \return a tensor viewing the computed values if they are stored in memory
     * owned by another object, a tensor owning a copy of these values otherwise
     */
    constexpr auto apply(const evaluator_result<EvaluatorType>&) const;
  };

  /*!
   * \brief perform consistency checks
   * \tparam TensorType: type of the tensorial object to be exposed
   * \tparam EvaluatorType: type of the underlying evaluator
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: tensor modifier
   * \return true if the values computed by the modifier are valid
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const TensorModifier<TensorType, EvaluatorType>&);
  /*!
   * \brief return the number of components
   * \tparam TensorType: type of the tensorial object to be exposed
   * \tparam EvaluatorType: type of the underlying evaluator
   * \param[in] e: tensor modifier
   * \return the number of components of the modifier, which is the number of
   * components of the tensorial type
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorModifier<TensorType, EvaluatorType>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX */
