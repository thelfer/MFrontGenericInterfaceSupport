/*!
 * \file   MGIS/Function/FixedSizeView.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX

#include "MGIS/Function/Space.hxx"
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
  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable = false>
  requires(N > 0) struct FixedSizeView {
    //! \brief value returned by non-const call operator
    using mutable_value_type =
        std::conditional_t<N == 1, real&, std::span<real, N>>;
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] values: function
     */
    static bool checkPreconditions(
        const FunctionView<Space, {}, is_mutable>&) noexcept;
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    FixedSizeView(const FunctionView<Space, {}, is_mutable>&);
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
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    mutable_value_type operator()(const element_index<Space>&) requires(
        is_mutable&& ElementSpaceConcept<Space> &&
        !(hasElementWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    mutable_value_type operator()(
        const element_workspace<Space>&,
        const element_index<
            Space>&) requires(is_mutable&& ElementSpaceConcept<Space>&&
                                  hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    mutable_value_type operator()(
        const cell_index<Space>,
        const quadrature_point_index<
            Space>) requires(is_mutable&& QuadratureSpaceConcept<Space> &&
                             (!hasCellWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    mutable_value_type operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>,
        const quadrature_point_index<
            Space>) requires(is_mutable&& QuadratureSpaceConcept<Space>&&
                                 hasCellWorkspace<Space>);

   private:
    //! \brief underlying function
    std::conditional_t<is_mutable,
                       FunctionView<Space, {}, is_mutable>,
                       const FunctionView<Space, {}, is_mutable>>
        function;
  };  // end of FixedSizeView

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space, bool is_mutable>
  auto view(const FunctionView<Space, {}, is_mutable>&) requires(N > 0);

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const Function<Space, dynamic_extent>&)  //
      requires((N > 0) && (N != dynamic_extent));
  /*!
   * \brief convert a function to a mutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space>
  auto view(Function<Space, dynamic_extent>&)  //
      requires((N > 0) && (N != dynamic_extent));

}  // end of namespace mgis::function

#include "MGIS/Function/FixedSizeView.ixx"

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_HXX */
