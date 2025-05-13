/*!
 * \file   MGIS/Function/TensorView.hxx
 * \brief  This file declares the TensorView class and the
 * transform function, as well as the as_fsarray, as_tvector,
 * as_tmatrix, as_stensor and as_tensor modifiers.
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TENSORVIEW_HXX

#include <span>
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of an immutable
   * function view as a fixed size span or a scalar
   *
   * \tparam Space: functional space
   * \tparam TensorType: type of the tensor
   * \tparam is_mutable: boolean stating if the call operators can return
   * mutable values
   */
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable = false>
  struct TensorView {
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
    TensorView(const FunctionView<Space, {}, is_mutable>&);
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
    auto operator()(const element_index<Space>&) requires(
        is_mutable&& ElementSpaceConcept<Space> &&
        !(hasElementWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(
        const element_workspace<Space>&,
        const element_index<
            Space>&) requires(is_mutable&& ElementSpaceConcept<Space>&&
                                  hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(
        const cell_index<Space>,
        const quadrature_point_index<
            Space>) requires(is_mutable&& QuadratureSpaceConcept<Space> &&
                             (!hasCellWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(
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
  };  // end of TensorView

}  // end of namespace mgis::function

#include "MGIS/Function/Tensors/TensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORVIEW_HXX */
