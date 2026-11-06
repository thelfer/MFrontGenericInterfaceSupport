/*!
 * \file   MGIS/Function/TFEL/QuantityView.hxx
 * \brief  File which defines the `QuantityView` class, a view returning the
 * values of a function with one component as TFEL's quantities, i.e. scalar
 * values associated with a physical unit.
 * \author Thomas Helfer
 * \date   15/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_HXX

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#include "TFEL/Math/qt.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/TFEL/General.hxx"

namespace mgis::function {

  namespace internals {

    /*!
     * \brief partial specialization allowing functions to return a mutable
     * reference to a quantity, i.e. a `qt_ref` object of `TFEL/Math`
     *
     * \tparam UnitType: unit of the quantity
     */
    template <::tfel::math::unit::UnitConcept UnitType>
    struct FunctionResultTypeTraits<::tfel::math::qt_ref<UnitType, real>> {
      //! \brief the traits are specialized for this type
      static constexpr auto is_specialized = true;
    };

  }  // namespace internals

  /*!
   * \brief a view returning the values of a function as quantities
   *
   * A quantity is a scalar value associated with a physical unit. The view is
   * built from a function with exactly one component. Its call operators
   * return the value of the function at a given integration point:
   *
   * - a mutable quantity (a `qt_ref` object of `TFEL/Math`), referencing the
   *   value of the function, for the non constant call operators,
   * - a read-only quantity (a `const_qt_ref` object of `TFEL/Math`),
   *   referencing the value of the function, for the constant call operators.
   *
   * Only the call operators which match the way the underlying function can be
   * called are available.
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam UnitType: unit of the quantities returned by the view
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct QuantityView
      : private PreconditionsChecker<QuantityView<FunctionType, UnitType>> {
    //! \brief type of the space on which the underlying function is defined
    using Space = function_space<FunctionType>;
    /*!
     * \brief a simple alias used to workaround what seems to be a bug in
     * gcc 16.x: the underlying function is passed to the constructors by value
     * if it is a lightweight view and by reference otherwise
     */
    using ConstructorArgumentType =
        std::conditional_t<LightweightViewConcept<FunctionType>,
                           FunctionType,
                           FunctionType&>;
    /*!
     * \brief method checking that the preconditions of the constructor are met
     * \param[in] eh: error handler used to report the errors
     * \param[in] values: function used to build the view
     * \return true if the number of components of the function is equal to one
     */
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             const FunctionType&);
    /*!
     * \brief constructor, checking the preconditions
     * \param[in] values: function used to build the view
     */
    constexpr QuantityView(ConstructorArgumentType);
    /*!
     * \brief constructor
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function used to build the view
     */
    template <bool doPreconditionsCheck>
    constexpr QuantityView(const PreconditionsCheck<doPreconditionsCheck>&,
                           ConstructorArgumentType);
    /*!
     * \brief perform consistency checks
     * \param[in] eh: error handler used to report the errors
     * \return true if the number of components of the underlying function is
     * equal to one
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying function is defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \return the number of components of the values returned by the view,
     * which is one
     */
    [[nodiscard]] constexpr size_type getNumberOfComponents() const
        MGIS_NOEXCEPT;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point as a
     * read-only quantity
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
     * \return the value of the function at the given integration point as a
     * read-only quantity
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
     * \return the value of the function at the given integration point as a
     * read-only quantity
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
     * \return the value of the function at the given integration point as a
     * read-only quantity
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
     * \return the value of the function at the given integration point as a
     * mutable quantity, which may be modified
     * \note only available if the underlying function can be called with an
     * element index and if the values it returns can be mapped to mutable
     * memory
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator returning mutable values
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     * \return the value of the function at the given integration point as a
     * mutable quantity, which may be modified
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
     * \return the value of the function at the given integration point as a
     * mutable quantity, which may be modified
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
     * \return the value of the function at the given integration point as a
     * mutable quantity, which may be modified
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
  };  // end of QuantityView

  /*!
   * \brief partial specialisation stating that a `QuantityView` is its own
   * view: it only stores a view of the underlying function, so copying it is
   * cheap
   *
   * \tparam FunctionType: type of the underlying function
   * \tparam UnitType: unit of the quantities returned by the view
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct LightweightViewTraits<QuantityView<FunctionType, UnitType>>
      : std::true_type {};

  /*!
   * \brief return the underlying space
   * \tparam FunctionType: type of the underlying function
   * \tparam UnitType: unit of the quantities returned by the view
   * \param[in] v: quantity view
   * \return the space on which the view is defined
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const QuantityView<FunctionType, UnitType>&);
  /*!
   * \brief perform consistency checks
   * \tparam FunctionType: type of the underlying function
   * \tparam UnitType: unit of the quantities returned by the view
   * \param[in] eh: error handler used to report the errors
   * \param[in] v: quantity view
   * \return true if the values returned by the view are valid
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const QuantityView<FunctionType, UnitType>&);
  /*!
   * \brief return the number of components
   * \tparam FunctionType: type of the underlying function
   * \tparam UnitType: unit of the quantities returned by the view
   * \param[in] v: quantity view
   * \return the number of components of the view, which is one
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const QuantityView<FunctionType, UnitType>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/QuantityView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_HXX */
