/*!
 * \file   MGIS/Function/FixedSizeView.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_HXX

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of an immutable
   * function view as a fixed size span or a scalar
   *
   * \tparam Space: functional space
   * \tparam N: size of the returned value
   */
  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) struct FixedSizeView {
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] values: function
     */
    static bool checkPreconditions(
        const ImmutableFunctionView<Space, {}>&) noexcept;
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    FixedSizeView(const ImmutableFunctionView<Space, {}>&);
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
        requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_workspace<Space>&,
                    const element_index<Space>&) const
        requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_workspace<Space>&,
                    const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>);

   private:
    //! \brief underlying function
    const ImmutableFunctionView<Space, {}> function;
  };  // end of FixedSizeView

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const ImmutableFunctionView<Space, {}>&) requires(N > 0);

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const Function<Space, dynamic_extent>&)  //
      requires((N > 0) && (N != dynamic_extent));

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeView.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_HXX */
