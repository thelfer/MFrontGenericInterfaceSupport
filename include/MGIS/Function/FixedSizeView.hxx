/*!
 * \file   MGIS/Function/FixedSizeView.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX

#include "MGIS/Contract.hxx"
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
  template <FunctionConcept FunctionType, size_type N>
  requires(N > 0) struct FixedSizeView
      : private PreconditionsChecker<FixedSizeView<FunctionType, N>> {
    //
    using Space = function_space<FunctionType>;
    //! \brief value returned by non-const call operator
    using mutable_value_type =
        std::conditional_t<N == 1, real&, std::span<real, N>>;
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] eh: error handler
     * \param[in] values: function
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const FunctionType&);
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    constexpr FixedSizeView(FunctionType&);
    /*!
     * \brief constructor
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function
     */
    template <bool doPreconditionsCheck>
    constexpr FixedSizeView(const PreconditionsCheck<doPreconditionsCheck>&,
                            FunctionType&);
    //! \brief perform consistency checks
    constexpr bool check(AbstractErrorHandler&) const;
    //! \brief return the underlying  space
    constexpr decltype(auto) getSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_workspace<Space>&,
                              const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_index<Space>&,
                              const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_workspace<Space>&,
                              const cell_index<Space>&,
                              const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr mutable_value_type operator()(const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr mutable_value_type operator()(const element_workspace<Space>&,
                                            const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr mutable_value_type operator()(
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr mutable_value_type operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionResultTypeMappable<FunctionType>));

   private:
    //! \brief underlying function
    FunctionType& function;
  };  // end of FixedSizeView

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, typename FunctionType>
  constexpr auto view(FunctionType&) requires(
      (N > 0) && (N != dynamic_extent) &&              //
      (FunctionConcept<std::decay_t<FunctionType>>)&&  //
      (!std::is_rvalue_reference_v<FunctionType>));

  template <FunctionConcept FunctionType, size_type N>
  constexpr decltype(auto) getSpace(const FixedSizeView<FunctionType, N>&);

    //! \brief allocate internal workspace
  template <FunctionConcept FunctionType, size_type N>
  constexpr void allocateWorkspace(FixedSizeView<FunctionType, N>&) noexcept;
  //! \return the number of components
  template <FunctionConcept FunctionType, size_type N>
  constexpr size_type getNumberOfComponents(
      const FixedSizeView<FunctionType, N>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeView.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX */
