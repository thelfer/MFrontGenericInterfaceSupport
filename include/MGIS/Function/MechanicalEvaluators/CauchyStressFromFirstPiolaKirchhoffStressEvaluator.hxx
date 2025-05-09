/*!
 * \file   CauchyStressFromFirstPiolaKirchhoffStressEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX
#define LIB_MGIS_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/BinaryOperationEvaluatorBase.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the Cauchy stress stress using evaluators of
   * the deformation gradient and the first Piola-Kirchhoff stress
   *
   * \tparam N: space dimension
   * \tparam DeformationGradientEvaluatorType: evaluator of the deformation
   * gradient
   * \tparam PK1EvaluatorType: evaluator of the first Piola-Kirchhoff stress
   */
  template <unsigned short N,
            EvaluatorConcept PK1EvaluatorType,
            EvaluatorConcept DeformationGradientEvaluatorType>
  requires((N == 1) || (N == 2) || (N == 3))  //
      struct CauchyStressFromFirstPiolaKirchhoffStressEvaluator
      : public BinaryOperationEvaluatorBase<
            CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
                N,
                PK1EvaluatorType,
                DeformationGradientEvaluatorType>,
            PK1EvaluatorType,
            DeformationGradientEvaluatorType> {
    //! \brief inheriting constructors
    using BinaryOperationEvaluatorBase<
        CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
            N,
            PK1EvaluatorType,
            DeformationGradientEvaluatorType>,
        PK1EvaluatorType,
        DeformationGradientEvaluatorType>::BinaryOperationEvaluatorBase;
    //! \brief a simple alias
    using Space = std::decay_t<
        decltype(std::declval<std::decay_t<PK1EvaluatorType>>().getSpace())>;
    //! \brief perform consistency checks
    bool check(Context&) const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * do the conversion
     * \param[in] pk1_values
     * \param[in] F_values
     */
    auto apply(const auto&, const auto&) const;
  };

  namespace internals {

    template <unsigned short N,
              EvaluatorConcept DeformationGradientEvaluatorType>
    requires((N == 1) || (N == 2) ||
             (N == 3)) struct from_pk1_to_cauchy_modifier {
      //
      from_pk1_to_cauchy_modifier(DeformationGradientEvaluatorType&&);
      //
      template <typename PK1EvaluatorType>
      auto operator()(PK1EvaluatorType&&) const
          requires(EvaluatorConcept<std::decay_t<PK1EvaluatorType>>);

     private:
      DeformationGradientEvaluatorType F;
    };

  }  // namespace internals

  template <unsigned short N, typename DeformationGradientEvaluatorType>
  auto from_pk1_to_cauchy(DeformationGradientEvaluatorType&&) requires(
      EvaluatorConcept<std::decay_t<DeformationGradientEvaluatorType>>);

}  // end of namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.ixx"

#endif /* LIB_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX */
