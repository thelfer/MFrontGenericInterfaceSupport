/*!
 * \file   MGIS/Function/FixedSizeModifier.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX

#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of a
   * function view as a fixed size span or a scalar
   *
   * \tparam Space: functional space
   * \tparam N: size of the returned value
   */
  template <EvaluatorConcept EvaluatorType, size_type N>
  requires(N > 0) struct FixedSizeModifier {
    //
    using Space = evaluator_space<EvaluatorType>;
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] e: evaluator
     */
    static bool checkPreconditions(const EvaluatorType&) noexcept;
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    FixedSizeModifier(const EvaluatorType&);
    //! \brief perform consistency checks
    bool check(Context&) const noexcept;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying  space
    const Space& getSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b1);
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_workspace<Space>&,
                    const element_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b2);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_index<Space>&,
                    const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b3);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_workspace<Space>&,
                    const cell_index<Space>&,
                    const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b4);

   private:
    //! \brief underlying function
    EvaluatorType evaluator;
  };  // end of FixedSizeModifier

  /*!
   * \brief create a fixed size view from an evaluator
   * \param[in] f: function
   */
  template <size_type N, typename EvaluatorType>
  auto view(const EvaluatorType&) requires(
      (N > 0) && (N != dynamic_extent) &&  //
      (EvaluatorConcept<std::decay_t<EvaluatorType>>));

  template <EvaluatorConcept EvaluatorType, size_type N>
  const auto& getSpace(const FixedSizeModifier<EvaluatorType, N>&);

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEMODIFIER_HXX */
