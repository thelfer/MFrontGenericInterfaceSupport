/*!
 * \file   MGIS/Function/TFEL/Tensors.ixx
 * \brief  Implementation of the tensor modifiers and of the customization
 * points used to reduce tensorial values computed by evaluators, declared in
 * `MGIS/Function/TFEL/Tensors.hxx`.
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

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORS_IXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORS_IXX

namespace mgis::function::customization_points {

  /*!
   * \brief partial specialization computing the absolute value of a tensorial
   * object component by component
   *
   * The reduction is delegated to the `std::span<const real, N>`
   * specialization, applied to the values stored contiguously by the tensor.
   * The returned result owns a copy of the reduced values.
   *
   * \tparam TensorType: type of the tensorial object, which must match the
   * `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct AbsoluteValue<TensorType> {
    /*!
     * \brief compute the absolute value of every component of a tensorial
     * object
     * \param[in] v: tensorial object considered
     * \return a tensorial object holding the absolute values of the components
     * of `v`
     */
    static constexpr auto exe(const TensorType& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return AbsoluteValue<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief partial specialization computing the absolute value of a view on a
   * tensorial object component by component
   *
   * \tparam TensorType: type of the tensorial object viewed, which must match
   * the `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct AbsoluteValue<tfel::math::View<const TensorType>> {
    /*!
     * \brief compute the absolute value of every component of a viewed
     * tensorial object
     * \param[in] v: view on the tensorial object considered
     * \return a tensorial object holding the absolute values of the components
     * of `v`
     */
    static constexpr auto exe(
        const tfel::math::View<const TensorType>& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return AbsoluteValue<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief partial specialization returning the largest component of a
   * tensorial object
   *
   * \tparam TensorType: type of the tensorial object, which must match the
   * `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct MaximumComponent<TensorType> {
    /*!
     * \brief return the largest component of a tensorial object
     * \param[in] v: tensorial object considered
     * \return the largest of its components
     */
    static constexpr real exe(const TensorType& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return MaximumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief partial specialization returning the largest component of a view on
   * a tensorial object
   *
   * \tparam TensorType: type of the tensorial object viewed, which must match
   * the `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct MaximumComponent<tfel::math::View<const TensorType>> {
    /*!
     * \brief return the largest component of a viewed tensorial object
     * \param[in] v: view on the tensorial object considered
     * \return the largest of its components
     */
    static constexpr real exe(
        const tfel::math::View<const TensorType>& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return MaximumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief partial specialization returning the smallest component of a
   * tensorial object
   *
   * \tparam TensorType: type of the tensorial object, which must match the
   * `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct MinimumComponent<TensorType> {
    /*!
     * \brief return the smallest component of a tensorial object
     * \param[in] v: tensorial object considered
     * \return the smallest of its components
     */
    static constexpr real exe(const TensorType& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return MinimumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief partial specialization returning the smallest component of a view on
   * a tensorial object
   *
   * \tparam TensorType: type of the tensorial object viewed, which must match
   * the `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct MinimumComponent<tfel::math::View<const TensorType>> {
    /*!
     * \brief return the smallest component of a viewed tensorial object
     * \param[in] v: view on the tensorial object considered
     * \return the smallest of its components
     */
    static constexpr real exe(
        const tfel::math::View<const TensorType>& v) MGIS_NOEXCEPT {
      constexpr auto N = compile_time_size<TensorType>;
      return MinimumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

}  // end of namespace mgis::function::customization_points

namespace mgis::function::internals {

  /*!
   * \brief expose the values of a function argument as a tensorial object
   *
   * \tparam TensorType: type of the tensorial object to be exposed
   * \tparam QualifiedFunctionArgumentType: type of the function argument, which
   * may be a qualified type
   * \param[in] f: function argument whose values are exposed
   * \return a `TensorView` viewing the values of the function argument
   */
  template <TensorConcept TensorType>
  template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType>
  constexpr auto tensor_modifier<TensorType>::operator()(
      QualifiedFunctionArgumentType&& f) const
      requires(internals::checkNumberOfComponentsCompatibility<
               QualifiedFunctionArgumentType,
               compile_time_size<TensorType>>()) {
    return TensorView<std::decay_t<QualifiedFunctionArgumentType>, TensorType>(
        std::forward<QualifiedFunctionArgumentType>(f));
  }  // end of operator()

  /*!
   * \brief expose the values computed by an evaluator as a tensorial object
   *
   * \tparam TensorType: type of the tensorial object to be exposed
   * \tparam EvaluatorType: type of the underlying evaluator
   * \param[in] e: evaluator whose values are exposed
   * \return a `TensorModifier` returning the values computed by the evaluator
   * as objects of type `TensorType`
   */
  template <TensorConcept TensorType>
  template <EvaluatorConcept EvaluatorType>
  constexpr auto tensor_modifier<TensorType>::operator()(const EvaluatorType& e)
      const requires((areTensorModifierRequirementsSatisfied<
                      TensorType,
                      std::decay_t<EvaluatorType>>)) {
    using Modifier = TensorModifier<TensorType, std::decay_t<EvaluatorType>>;
    return Modifier(e);
  }  // end of operator()

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief expose the values of a function argument as a tensorial object, the
   * definition of the `operator|` used by expressions such as `f |
   * as_stensor<3>`
   *
   * \tparam QualifiedFunctionArgumentType: type of the function argument, which
   * may be a qualified type
   * \tparam TensorType: type of the tensorial object to be exposed
   * \param[in] f: function argument whose values are exposed
   * \param[in] m: modifier stating the type of the tensorial object to be
   * exposed
   * \return a `TensorView` viewing the values of the function argument
   */
  template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType,
            TensorConcept TensorType>
  constexpr auto operator|(QualifiedFunctionArgumentType&& f,
                           const internals::tensor_modifier<TensorType>& m)  //
      requires(internals::checkNumberOfComponentsCompatibility<
               QualifiedFunctionArgumentType,
               compile_time_size<TensorType>>()) {
    return m(std::forward<QualifiedFunctionArgumentType>(f));
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORS_IXX */
