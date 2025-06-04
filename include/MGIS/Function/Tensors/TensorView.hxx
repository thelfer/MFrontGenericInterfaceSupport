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
#include "MGIS/Contract.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of an immutable
   * function view as a fixed size span or a scalar
   *
   * \tparam FunctionType: underlying function type
   * \tparam TensorType: type of the tensor
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  struct TensorView
      : private PreconditionsChecker<TensorView<FunctionType, TensorType>> {
    //
    static_assert(number_of_components<FunctionType> == dynamic_extent
                      ? true
                      : compile_time_size<TensorType> ==
                            number_of_components<FunctionType>);
    //
    using Space = function_space<FunctionType>;
    /*!
     * \brief method checking that the precondition of the constructor are met.
     * \param[in] values: function
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const FunctionType&) noexcept;
    /*!
     * \brief constructor
     * \param[in] values: function
     */
    constexpr TensorView(FunctionType&);
    /*!
     * \brief constructor
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function
     */
    template <bool doPreconditionsCheck>
    constexpr TensorView(const PreconditionsCheck<doPreconditionsCheck>&,
                         FunctionType&);
    //! \brief perform consistency checks
    bool check(Context&) const noexcept;
    //! \brief allocate internal workspace
    constexpr void allocateWorkspace();
    //! \brief return the underlying  space
    constexpr const Space& getSpace() const;
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
    constexpr auto operator()(const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_workspace<Space>&,
                              const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_index<Space>&,
                              const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_workspace<Space>&,
                              const cell_index<Space>&,
                              const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionResultTypeMappable<FunctionType>));

   private:
    //! \brief underlying function
    FunctionType& function;
  };  // end of TensorView

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr decltype(auto) getSpace(
      const TensorView<FunctionType, TensorType>&);

}  // end of namespace mgis::function

#include "MGIS/Function/Tensors/TensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORVIEW_HXX */
