/*!
 * \file   MGIS/Function/TFEL/TensorView.hxx
 * \brief  This file declares the TensorView class and the
 * transform function, as well as the as_fsarray, as_tvector,
 * as_tmatrix, as_stensor and as_tensor modifiers.
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
#include "MGIS/Function/TFEL/TensorConcept.hxx"

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
                                             const FunctionType&);
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
    constexpr bool check(AbstractErrorHandler&) const;
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
    function_view<FunctionType> function;
  };  // end of TensorView

  //! \brief partial specialisation
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  struct LightweightViewTraits<TensorView<FunctionType, TensorType>>
      : std::true_type {};

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const TensorView<FunctionType, TensorType>&);
  //! \brief perform consistency checks
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&, const TensorView<FunctionType, TensorType>&);
  //! \return the number of components
  template <FunctionConcept FunctionType, TensorConcept TensorType>
  [[nodiscard]] constexpr size_type getNumberOfComponents(
      const TensorView<FunctionType, TensorType>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorView.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_HXX */
