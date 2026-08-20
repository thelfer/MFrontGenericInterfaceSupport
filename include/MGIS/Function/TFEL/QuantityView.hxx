/*!
 * \file   MGIS/Function/TFEL/QuantityView.hxx
 * \brief
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

namespace mgis::function {

  namespace internals {

    //! \brief partial specialization for reference to a std::array
    template <::tfel::math::unit::UnitConcept UnitType>
    struct FunctionResultTypeTraits<::tfel::math::qt_ref<UnitType, real>> {
      static constexpr auto is_specialized = true;
    };

  }  // namespace internals

  /*!
   * \brief a modifier returning the values of a
   * function view as a quantity
   *
   * \tparam Space: functional space
   */
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct QuantityView
      : private PreconditionsChecker<QuantityView<FunctionType, UnitType>> {
    //
    using Space = function_space<FunctionType>;
    //! \brief a simple alias used to workaround what seems to be a bug in gcc 16.x
    using ConstructorArgumentType =
        std::conditional_t<LightweightViewConcept<FunctionType>,
                           FunctionType,
                           FunctionType&>;
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
    constexpr QuantityView(ConstructorArgumentType);
    /*!
     * \brief constructor
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] values: function
     */
    template <bool doPreconditionsCheck>
    constexpr QuantityView(const PreconditionsCheck<doPreconditionsCheck>&,
                           ConstructorArgumentType);
    //! \brief perform consistency checks
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    //! \brief return the underlying  space
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    //! \return the number of components
    [[nodiscard]] constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] wk: element workspace
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionConstResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b1) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b2) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b3) &&
                 (isFunctionResultTypeMappable<FunctionType>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&)  //
        requires((internals::FunctionResultQuery<FunctionType>::b4) &&
                 (isFunctionResultTypeMappable<FunctionType>));

   private:
    //! \brief underlying view
    function_view<FunctionType> function;
  };  // end of QuantityView

  //! \brief partial specialisation
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  struct LightweightViewTraits<QuantityView<FunctionType, UnitType>>
      : std::true_type {};

  //! \return the underlying space
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const QuantityView<FunctionType, UnitType>&);
  //! \brief perform consistency checks
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const QuantityView<FunctionType, UnitType>&);
  //! \return the number of components
  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const QuantityView<FunctionType, UnitType>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/QuantityView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_HXX */
