/*!
 * \file   MGIS/Function/Mechanics.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_MECHANICS_IXX
#define LIB_MGIS_FUNCTION_MECHANICS_IXX

#ifdef MGIS_HAVE_TFEL

namespace mgis::function::internals {

  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  requires(
      std::is_convertible_v<evaluator_result<StiffnessEvaluator>,
                            tfel::material::tangent_operator<ResultFlag,
                                                             N,
                                                             ::mgis::real>>)  //
      struct ConvertFiniteStrainStiffnessEvaluator {
    //! \brief a simple alias
    using Space = evaluator_space<StiffnessEvaluator>;
    /*!
     * \brief constructor
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
    //! \brief perform consistency checks
    [[nodiscard]] constexpr bool check(AbstractErrorHandler& ctx) const {
      return ((checkMatchingSpaces(ctx, this->K, this->F0)) &&
              (checkMatchingSpaces(ctx, this->K, this->F1)) &&
              (checkMatchingSpaces(ctx, this->K, this->s)) &&
              (internals::disambiguateCheck(ctx, this->K)) &&
              (internals::disambiguateCheck(ctx, this->F0)) &&
              (internals::disambiguateCheck(ctx, this->F1)) &&
              (internals::disambiguateCheck(ctx, this->s)));
    }

    //! \brief allocate internal workspace
    constexpr void allocateWorkspace() {
      internals::disambiguateAllocateWorkspace(this->K);
      internals::disambiguateAllocateWorkspace(this->F0);
      internals::disambiguateAllocateWorkspace(this->F1);
      internals::disambiguateAllocateWorkspace(this->s);
    }
    //! \brief return the underlying space
    [[nodiscard]] constexpr decltype(auto) getSpace() const {
      return internals::disambiguateGetSpace(this->K);
    }
    /*!
     * \brief call operator
     * \param[in] i: element index
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
     * \brief call operator
     * \param[in] wk: element index
     * \param[in] i: element index
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
          this->K(wk, i), this->F0(wk, i), this->F1(wk, i), this->s(wk, i));
    }
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
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
          this->K(e, i), this->F0(e, i), this->F1(e, i), this->s(e, i));
    }
    /*!
     * \brief call operator
     * \param[in] wk: workspace
     * \param[in] e: cell index
     * \param[in] i: integration point index
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
          this->K(wk, e, i), this->F0(wk, e, i), this->F1(wk, e, i),
          this->s(wk, e, i));
    }

   private:
    StiffnessEvaluator K;
    DeformationGradientEvaluatorType0 F0;
    DeformationGradientEvaluatorType1 F1;
    CauchyStressEvaluatorType s;
  };

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

  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            FourthOrderTensorEvaluatorConcept StiffnessEvaluator,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  constexpr void allocateWorkspace(
      ConvertFiniteStrainStiffnessEvaluator<N,
                                            ResultFlag,
                                            SourceFlag,
                                            StiffnessEvaluator,
                                            DeformationGradientEvaluatorType0,
                                            DeformationGradientEvaluatorType1,
                                            CauchyStressEvaluatorType>& e) {
    return e.allocateWorkspace();
  }  // end of allocateWorkspace

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
          CauchyStressEvaluatorType>&) noexcept {
    const auto value =
        tfel::material::tangent_operator<SourceFlag, N, ::mgis::real>{};
    return value.size();
  }  // end of getNumberOfComments

  template <unsigned short N,
            FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  struct ConvertFiniteStrainStiffnessModifierGenerator {
    constexpr ConvertFiniteStrainStiffnessModifierGenerator(
        const DeformationGradientEvaluatorType0& F0_,
        const DeformationGradientEvaluatorType1& F1_,
        const CauchyStressEvaluatorType& s_)
        : F0(F0_),
          F1(F1_),
          s(s_) {} /* end of ConvertFiniteStrainStiffnessModifierGenerator */

    using Tag = ::mgis::function::EvaluatorModifierTag;

    template <FourthOrderTensorEvaluatorConcept StiffnessEvaluator>
    [[nodiscard]] constexpr auto
    operator()(const StiffnessEvaluator& K) const requires(
        std::is_convertible_v<
            evaluator_result<StiffnessEvaluator>,
            tfel::material::tangent_operator<ResultFlag, N, ::mgis::real>>) {
      return ConvertFiniteStrainStiffnessEvaluator<
          N, ResultFlag, SourceFlag, StiffnessEvaluator,
          DeformationGradientEvaluatorType0, DeformationGradientEvaluatorType1,
          CauchyStressEvaluatorType>(K, F0, F1, s);
    }  // end of operator()

   private:
    const DeformationGradientEvaluatorType0 F0;
    const DeformationGradientEvaluatorType1 F1;
    const CauchyStressEvaluatorType s;
  };

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

#endif MGIS_HAVE_TFEL

#endif /* LIB_MGIS_FUNCTION_MECHANICS_IXX */
