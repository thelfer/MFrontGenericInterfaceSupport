/*!
 * \file   MGIS/Function/TFEL/Mechanics.ixx
 * \brief  Implementation of the finite strain tangent operator conversion
 * declared in `MGIS/Function/TFEL/Mechanics.hxx`: the
 * `ConvertFiniteStrainStiffnessEvaluator` evaluator, the associated modifier
 * generator and the `convert_finite_strain_stiffness` function.
 * \author Thomas Helfer
 * \date   02/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_MECHANICS_IXX
#define LIB_MGIS_FUNCTION_TFEL_MECHANICS_IXX

namespace mgis::function::internals {

  /*!
   * \brief evaluator computing a finite strain tangent operator from a tangent
   * operator of another kind, given two deformation gradients and the Cauchy
   * stress
   *
   * The conversion is performed by
   * `tfel::material::convert<ResultFlag, SourceFlag, N>`, applied to the values
   * returned by the four underlying evaluators at a given integration point.
   * Only the call operators which match the way all the underlying evaluators
   * can be called are available.
   *
   * \tparam N: space dimension (1, 2 or 3)
   * \tparam ResultFlag: kind of the tangent operator returned by the evaluator
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam StiffnessEvaluator: type of the evaluator providing the tangent
   * operator to convert
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \note only available if the values returned by the stiffness evaluator are
   * convertible to a `tfel::material::tangent_operator` of kind `SourceFlag`
   */
  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  requires(
      std::is_convertible_v<evaluator_result<StiffnessEvaluator>,
                            tfel::material::tangent_operator<SourceFlag,
                                                             N,
                                                             ::mgis::real>>)  //
      struct ConvertFiniteStrainStiffnessEvaluator {
    //! \brief space on which the evaluator is defined
    using Space = evaluator_space<StiffnessEvaluator>;
    /*!
     * \brief constructor
     * \param[in] K_: evaluator providing the tangent operator to convert
     * \param[in] F0_: evaluator providing the first deformation gradient
     * \param[in] F1_: evaluator providing the second deformation gradient
     * \param[in] s_: evaluator providing the Cauchy stress
     */
    constexpr ConvertFiniteStrainStiffnessEvaluator(
        const StiffnessEvaluator& K_,
        const DeformationGradientEvaluatorType0& F0_,
        const DeformationGradientEvaluatorType1& F1_,
        const CauchyStressEvaluatorType& s_)
        : K(K_),
          F0(F0_),
          F1(F1_),
          s(s_) {}  // end of ConvertFiniteStrainStiffnessModifierGenerator
    //! \brief copy constructor
    constexpr ConvertFiniteStrainStiffnessEvaluator(
        const ConvertFiniteStrainStiffnessEvaluator&) = default;
    //! \brief move constructor
    constexpr ConvertFiniteStrainStiffnessEvaluator(
        ConvertFiniteStrainStiffnessEvaluator&&) = default;
    /*!
     * \brief perform consistency checks
     * \param[in] ctx: error handler used to report the errors
     * \return true if the four underlying evaluators are valid and are defined
     * on matching spaces
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler& ctx) const {
      return ((checkMatchingSpaces(ctx, this->K, this->F0)) &&
              (checkMatchingSpaces(ctx, this->K, this->F1)) &&
              (checkMatchingSpaces(ctx, this->K, this->s)) &&
              (internals::disambiguateCheck(ctx, this->K)) &&
              (internals::disambiguateCheck(ctx, this->F0)) &&
              (internals::disambiguateCheck(ctx, this->F1)) &&
              (internals::disambiguateCheck(ctx, this->s)));
    }

    /*!
     * \brief return the underlying space
     * \return the space of the stiffness evaluator
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const {
      return internals::disambiguateGetSpace(this->K);
    }
    /*!
     * \brief call operator computing the converted tangent operator at a given
     * element
     * \param[in] i: element index
     * \return the tangent operator of kind `ResultFlag` computed from the
     * values returned by the underlying evaluators
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>& i) const
        requires(
            (internals::EvaluatorResultQuery<StiffnessEvaluator>::b1) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType0>::b1) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType1>::b1) &&
            (internals::EvaluatorResultQuery<CauchyStressEvaluatorType>::b1)) {
      return tfel::material::convert<ResultFlag, SourceFlag, N>(
          this->K(i), tfel::math::tensor<N, ::mgis::real>(this->F0(i)),
          tfel::math::tensor<N, ::mgis::real>(this->F1(i)),
          tfel::math::stensor<N, ::mgis::real>(this->s(i)));
    }
    /*!
     * \brief call operator computing the converted tangent operator at a given
     * element, using an element workspace
     * \param[in] wk: workspace associated with the element
     * \param[in] i: element index
     * \return the tangent operator of kind `ResultFlag` computed from the
     * values returned by the underlying evaluators
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>& wk,
                                            const element_index<Space>& i) const
        requires(
            (internals::EvaluatorResultQuery<StiffnessEvaluator>::b2) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType0>::b2) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType1>::b2) &&
            (internals::EvaluatorResultQuery<CauchyStressEvaluatorType>::b2)) {
      return tfel::material::convert<ResultFlag, SourceFlag, N>(
          this->K(wk, i), tfel::math::tensor<N, ::mgis::real>(this->F0(wk, i)),
          tfel::math::tensor<N, ::mgis::real>(this->F1(wk, i)),
          tfel::math::stensor<N, ::mgis::real>(this->s(wk, i)));
    }
    /*!
     * \brief call operator computing the converted tangent operator at a given
     * quadrature point
     * \param[in] e: cell index
     * \param[in] i: integration point index
     * \return the tangent operator of kind `ResultFlag` computed from the
     * values returned by the underlying evaluators
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>& e,
        const quadrature_point_index<Space>& i) const
        requires(
            (internals::EvaluatorResultQuery<StiffnessEvaluator>::b3) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType0>::b3) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType1>::b3) &&
            (internals::EvaluatorResultQuery<CauchyStressEvaluatorType>::b3)) {
      return tfel::material::convert<ResultFlag, SourceFlag, N>(
          this->K(e, i), tfel::math::tensor<N, ::mgis::real>(this->F0(e, i)),
          tfel::math::tensor<N, ::mgis::real>(this->F1(e, i)),
          tfel::math::stensor<N, ::mgis::real>(this->s(e, i)));
    }
    /*!
     * \brief call operator computing the converted tangent operator at a given
     * quadrature point, using a cell workspace
     * \param[in] wk: workspace associated with the cell
     * \param[in] e: cell index
     * \param[in] i: integration point index
     * \return the tangent operator of kind `ResultFlag` computed from the
     * values returned by the underlying evaluators
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>& wk,
        const cell_index<Space>& e,
        const quadrature_point_index<Space>& i) const
        requires(
            (internals::EvaluatorResultQuery<StiffnessEvaluator>::b4) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType0>::b4) &&
            (internals::EvaluatorResultQuery<
                DeformationGradientEvaluatorType1>::b4) &&
            (internals::EvaluatorResultQuery<CauchyStressEvaluatorType>::b4)) {
      return tfel::material::convert<ResultFlag, SourceFlag, N>(
          this->K(wk, e, i),
          tfel::math::tensor<N, ::mgis::real>(this->F0(wk, e, i)),
          tfel::math::tensor<N, ::mgis::real>(this->F1(wk, e, i)),
          tfel::math::stensor<N, ::mgis::real>(this->s(wk, e, i)));
    }

   private:
    //! \brief evaluator providing the tangent operator to convert
    StiffnessEvaluator K;
    //! \brief evaluator providing the first deformation gradient
    DeformationGradientEvaluatorType0 F0;
    //! \brief evaluator providing the second deformation gradient
    DeformationGradientEvaluatorType1 F1;
    //! \brief evaluator providing the Cauchy stress
    CauchyStressEvaluatorType s;
  };

  /*!
   * \brief perform consistency checks on a finite strain stiffness conversion
   * evaluator
   *
   * \tparam N: space dimension (1, 2 or 3)
   * \tparam ResultFlag: kind of the tangent operator returned by the evaluator
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam StiffnessEvaluator: type of the evaluator providing the tangent
   * operator to convert
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: conversion evaluator
   * \return true if the conversion evaluator is valid
   */
  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler& eh,
      const ConvertFiniteStrainStiffnessEvaluator<
          N,
          ResultFlag,
          SourceFlag,
          StiffnessEvaluator,
          DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1,
          CauchyStressEvaluatorType>& e) {
    return e.check(eh);
  }  // end of check

  /*!
   * \brief return the number of components of a finite strain stiffness
   * conversion evaluator
   *
   * \tparam N: space dimension (1, 2 or 3)
   * \tparam ResultFlag: kind of the tangent operator returned by the evaluator
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam StiffnessEvaluator: type of the evaluator providing the tangent
   * operator to convert
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \param[in] e: conversion evaluator
   * \return the number of components of a tangent operator of kind `ResultFlag`
   * and dimension `N`
   */
  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  [[nodiscard]] constexpr mgis::size_type getNumberOfComponents(
      const ConvertFiniteStrainStiffnessEvaluator<
          N,
          ResultFlag,
          SourceFlag,
          StiffnessEvaluator,
          DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1,
          CauchyStressEvaluatorType>&) MGIS_NOEXCEPT {
    const auto value =
        tfel::material::tangent_operator<ResultFlag, N, ::mgis::real>{};
    return value.size();
  }  // end of getNumberOfComments

  /*!
   * \brief modifier which, applied to an evaluator returning a tangent operator
   * of kind `SourceFlag`, builds a `ConvertFiniteStrainStiffnessEvaluator`
   * converting it into a tangent operator of kind `ResultFlag`
   *
   * \tparam N: space dimension (1, 2 or 3)
   * \tparam ResultFlag: kind of the tangent operator to be returned
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   */
  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  struct ConvertFiniteStrainStiffnessModifierGenerator {
    /*!
     * \brief constructor
     * \param[in] F0_: evaluator providing the first deformation gradient
     * \param[in] F1_: evaluator providing the second deformation gradient
     * \param[in] s_: evaluator providing the Cauchy stress
     */
    constexpr ConvertFiniteStrainStiffnessModifierGenerator(
        const DeformationGradientEvaluatorType0& F0_,
        const DeformationGradientEvaluatorType1& F1_,
        const CauchyStressEvaluatorType& s_)
        : F0(F0_),
          F1(F1_),
          s(s_) {} /* end of ConvertFiniteStrainStiffnessModifierGenerator */

    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;

    /*!
     * \brief build the evaluator converting a tangent operator
     * \tparam StiffnessEvaluator: type of the evaluator providing the tangent
     * operator to convert
     * \param[in] K: evaluator providing the tangent operator to convert
     * \return a `ConvertFiniteStrainStiffnessEvaluator`
     * \note only available if the values returned by `K` are convertible to a
     * `tfel::material::tangent_operator` of kind `SourceFlag`
     */
    template <FourthOrderTensorEvaluatorConcept StiffnessEvaluator>
    [[nodiscard]] constexpr auto
    operator()(const StiffnessEvaluator& K) const requires(
        std::is_convertible_v<
            evaluator_result<StiffnessEvaluator>,
            tfel::material::tangent_operator<SourceFlag, N, ::mgis::real>>) {
      return ConvertFiniteStrainStiffnessEvaluator<
          N, ResultFlag, SourceFlag, StiffnessEvaluator,
          DeformationGradientEvaluatorType0, DeformationGradientEvaluatorType1,
          CauchyStressEvaluatorType>(K, F0, F1, s);
    }  // end of operator()

   private:
    //! \brief evaluator providing the first deformation gradient
    const DeformationGradientEvaluatorType0 F0;
    //! \brief evaluator providing the second deformation gradient
    const DeformationGradientEvaluatorType1 F1;
    //! \brief evaluator providing the Cauchy stress
    const CauchyStressEvaluatorType s;
  };

  /*!
   * \brief return the space of a finite strain stiffness conversion evaluator
   *
   * \tparam N: space dimension (1, 2 or 3)
   * \tparam ResultFlag: kind of the tangent operator returned by the evaluator
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam StiffnessEvaluator: type of the evaluator providing the tangent
   * operator to convert
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \param[in] e: conversion evaluator
   * \return the space on which the evaluator is defined
   */
  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  [[nodiscard]] decltype(auto) getSpace(
      const ConvertFiniteStrainStiffnessEvaluator<
          N,
          ResultFlag,
          SourceFlag,
          StiffnessEvaluator,
          DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1,
          CauchyStressEvaluatorType>& e) {
    return e.getSpace();
  }  // end of getSpace

}  // namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief definition of the `convert_finite_strain_stiffness` function
   * declared in `MGIS/Function/TFEL/Mechanics.hxx`
   *
   * The space dimension `N` is deduced from that of the first deformation
   * gradient evaluator; the second deformation gradient and the Cauchy stress
   * must have the same dimension, which is enforced by `static_assert`.
   *
   * \tparam ResultFlag: kind of the tangent operator to be returned
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \param[in] F0: evaluator used to compute the first deformation gradient
   * \param[in] F1: evaluator used to compute the second deformation gradient
   * \param[in] s: evaluator used to compute the Cauchy stress
   * \return a modifier converting the tangent operators it is applied to
   */
  template <FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  constexpr auto convert_finite_strain_stiffness(
      const DeformationGradientEvaluatorType0& F0,
      const DeformationGradientEvaluatorType1& F1,
      const CauchyStressEvaluatorType& s) {
    if constexpr (Tensor1DEvaluatorConcept<DeformationGradientEvaluatorType0>) {
      static_assert(
          Tensor1DEvaluatorConcept<DeformationGradientEvaluatorType1>);
      static_assert(Stensor1DEvaluatorConcept<CauchyStressEvaluatorType>);
      return internals::ConvertFiniteStrainStiffnessModifierGenerator<
          1u, ResultFlag, SourceFlag, DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1, CauchyStressEvaluatorType>(F0, F1,
                                                                        s);
    } else if constexpr (Tensor2DEvaluatorConcept<
                             DeformationGradientEvaluatorType0>) {
      static_assert(
          Tensor2DEvaluatorConcept<DeformationGradientEvaluatorType1>);
      static_assert(Stensor2DEvaluatorConcept<CauchyStressEvaluatorType>);
      return internals::ConvertFiniteStrainStiffnessModifierGenerator<
          2u, ResultFlag, SourceFlag, DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1, CauchyStressEvaluatorType>(F0, F1,
                                                                        s);
    } else {
      static_assert(
          Tensor3DEvaluatorConcept<DeformationGradientEvaluatorType1>);
      static_assert(Stensor3DEvaluatorConcept<CauchyStressEvaluatorType>);
      return internals::ConvertFiniteStrainStiffnessModifierGenerator<
          3u, ResultFlag, SourceFlag, DeformationGradientEvaluatorType0,
          DeformationGradientEvaluatorType1, CauchyStressEvaluatorType>(F0, F1,
                                                                        s);
    }
  }  // end of convert_finite_strain_stiffness

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_MECHANICS_IXX */
