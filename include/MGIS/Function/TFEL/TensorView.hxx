/*!
 * \file   MGIS/Function/TFEL/TensorView.hxx
 * \brief  File which declares the `TensorView` class, an evaluator returning
 * the values of a function as TFEL's tensorial objects.
 * \author Thomas Helfer
 * \date   09/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_HXX

#include <span>
#include "MGIS/Contract.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the values of a function as TFEL's
   * tensorial objects
   *
   * The view is built from a function whose number of components must be
   * equal to the number of components of `TensorType`. Its call operators
   * return the values of the function at a given integration point, mapped on
   * the memory of the function using `tfel::math::map`:
   *
   * - a read-only tensor, for the constant call operators,
   * - a mutable tensor, which may be modified, for the non constant call
   *   operators.
   *
   * Only the call operators which match the way the underlying function can be
   * called are available.
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam TensorType: type of the tensorial objects returned by the view,
   * which must match the `TensorConcept` concept
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  struct TensorView
      : private PreconditionsChecker<TensorView<FunctionType, TensorType>> {
    //
    static_assert(number_of_components<FunctionType> == dynamic_extent
                      ? true
                      : compile_time_size<TensorType> ==
                            number_of_components<FunctionType>);
    //! \brief type of the space on which the underlying function is defined
    using Space = function_space<FunctionType>;
    /*!
     * \brief a simple alias used to workaround what seems to be a bug in
     * gcc 16.x: the underlying function is passed to the constructors by
     * value if it is a lightweight view and by reference otherwise
     */
    using ConstructorArgumentType =
        std::conditional_t<LightweightViewConcept<FunctionType>,
                           FunctionType,
                           FunctionType&>;
    /*!
     * \brief method checking that the preconditions of the constructor are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] values: function used to build the view
     * \return true if the number of components of the function is equal to the
     * number of components of the tensorial type
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const FunctionType&);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] values: function used to build the view
     */
    constexpr TensorView(ConstructorArgumentType);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function used to build the view
     */
    template <bool doPreconditionsCheck>
    constexpr TensorView(const PreconditionsCheck<doPreconditionsCheck>&,
                         ConstructorArgumentType);
    /*!
     * \brief perform consistency checks
     * \param[in] eh: error handler used to report the errors
     * \return true if the number of components of the underlying function is
     * equal to the number of components of the tensorial type
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying function is defined
     */
    [[nodiscard]] constexpr const Space& getSpace() const;
    /*!
     * \return the number of components of the values returned by the view,
     * which is the number of components of the tensorial type
     */
    [[nodiscard]] constexpr size_type getNumberOfComponents() const
        MGIS_NOEXCEPT;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     * \return the values of the function at the given integration point as a
     * read-only tensor
     * \note only available if the underlying function can be called with an
     * element index and if the values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the values of the function at the given integration point as a
     * read-only tensor
     * \note only available if the underlying function can be called with an
     * element workspace and an element index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values of the function at the given integration point as a
     * read-only tensor
     * \note only available if the underlying function can be called with a
     * cell index and a quadrature point index and if the values it returns can
     * be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values of the function at the given integration point as a
     * read-only tensor
     * \note only available if the underlying function can be called with a
     * cell workspace, a cell index and a quadrature point index and if the
     * values it returns can be mapped to memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] i: integration point index
     * \return the values of the function at the given integration point as a
     * mutable tensor, which may be modified
     * \note only available if the underlying function can be called with an
     * element index and if the values it returns can be mapped to mutable
     * memory
     */
    constexpr auto operator()(const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the values of the function at the given integration point as a
     * mutable tensor, which may be modified
     * \note only available if the underlying function can be called with an
     * element workspace and an element index and if the values it returns can
     * be mapped to mutable memory
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values of the function at the given integration point as a
     * mutable tensor, which may be modified
     * \note only available if the underlying function can be called with a
     * cell index and a quadrature point index and if the values it returns can
     * be mapped to mutable memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the values of the function at the given integration point as a
     * mutable tensor, which may be modified
     * \note only available if the underlying function can be called with a
     * cell workspace, a cell index and a quadrature point index and if the
     * values it returns can be mapped to mutable memory
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionResultTypeMappable<FunctionType>));

   private:
    //! \brief underlying view of the function
    function_view<FunctionType> function;
  };  // end of TensorView

  /*!
   * \brief partial specialisation stating that a `TensorView` is its own
   * view: it only stores a view of the underlying function, so copying it is
   * cheap
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam TensorType: type of the tensorial objects returned by the view
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  struct LightweightViewTraits<TensorView<FunctionType, TensorType>>
      : std::true_type {};

  /*!
   * \brief return the underlying space
   * \tparam FunctionType: type of the underlying function
   * \tparam TensorType: type of the tensorial objects returned by the view
   * \param[in] v: tensor view
   * \return the space on which the view is defined
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const TensorView<FunctionType, TensorType>&);
  /*!
   * \brief perform consistency checks
   * \tparam FunctionType: type of the underlying function
   * \tparam TensorType: type of the tensorial objects returned by the view
   * \param[in] eh: error handler used to report the errors
   * \param[in] v: tensor view
   * \return true if the values returned by the view are valid
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const TensorView<FunctionType, TensorType>&);
  /*!
   * \brief return the number of components
   * \tparam FunctionType: type of the underlying function
   * \tparam TensorType: type of the tensorial objects returned by the view
   * \param[in] v: tensor view
   * \return the number of components of the view, which is the number of
   * components of the tensorial type
   */
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const TensorView<FunctionType, TensorType>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_HXX */
