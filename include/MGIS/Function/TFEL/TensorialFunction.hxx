/*!
 * \file   MGIS/Function/TFEL/TensorialFunction.hxx
 * \brief  File which defines the `TensorialFunction` class, a function which
 * owns its values and returns them as TFEL's tensorial objects, and the
 * aliases for the most common tensorial functions.
 * \author Thomas Helfer
 * \date   28/08/2025
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

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_HXX

#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"

namespace mgis::function {

  /*!
   * \brief default implementation of a function whose values are TFEL's
   * tensorial objects
   *
   * A `TensorialFunction` stores its values in a `Function`, whose number of
   * components is the minimal size of the array of values required to store a
   * `TensorType`, as computed by the
   * `tfel::math::getUnderlyingArrayMinimalSize` function, and derives from the
   * corresponding `TensorView`, so it can be called like a function and returns
   * its values as tensorial objects viewing its storage.
   *
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function,
   * which must match the `TensorConcept` concept
   *
   * \note the data stride is equal to the data size
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      struct TensorialFunction
      : private PreconditionsChecker<TensorialFunction<Space, TensorType>>,
        private Function<Space,
                         tfel::math::getUnderlyingArrayMinimalSize<
                             typename TensorType::indexing_policy>()>,
        public TensorView<Function<Space,
                                   tfel::math::getUnderlyingArrayMinimalSize<
                                       typename TensorType::indexing_policy>()>,
                          TensorType> {
    //! \brief type of the function which owns the values
    using UnderlyingFunctionType =
        Function<Space,
                 tfel::math::getUnderlyingArrayMinimalSize<
                     typename TensorType::indexing_policy>()>;
    //! \brief type of the view through which the values are returned
    using UnderlyingTensorView =
        TensorView<Function<Space,
                            tfel::math::getUnderlyingArrayMinimalSize<
                                typename TensorType::indexing_policy>()>,
                   TensorType>;
    /*!
     * \brief check that the preconditions to build the function are met
     * \param[in] eh: unused error handler
     * \param[in] s: unused space
     * \return always true: no check is required
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&, const Space&);
    /*!
     * \brief constructor from a space, checking the preconditions
     * \param[in] s: space on which the function is defined
     */
    constexpr TensorialFunction(const Space&);
    /*!
     * \brief constructor from a space
     * \tparam doPreconditionsCheck: boolean stating if the preconditions must
     * be checked
     * \param[in] pcheck: object stating if preconditions must be checked
     * \param[in] s: space on which the function is defined
     */
    template <bool doPreconditionsCheck>
    constexpr TensorialFunction(const PreconditionsCheck<doPreconditionsCheck>&,
                                const Space&);
    /*!
     * \brief move constructor: moves the values of the function and rebinds the
     * view to them
     * \param[in] f: function to move
     */
    constexpr TensorialFunction(TensorialFunction&&);
    /*!
     * \brief copy constructor: copies the values of the function and rebinds
     * the view to them
     * \param[in] f: function to copy
     */
    constexpr TensorialFunction(const TensorialFunction&);
    /*!
     * \brief return a view of the function
     * \return a mutable view of the values of the function
     */
    constexpr UnderlyingTensorView& view() & MGIS_NOEXCEPT;
    /*!
     * \brief return a view of the function
     * \return a read-only view of the values of the function
     */
    constexpr const UnderlyingTensorView& view() const& MGIS_NOEXCEPT;
    //
    //! \brief creating a view of a temporary function is forbidden
    auto view() && = delete;
    //
    //! \brief inherit the `getNumberOfComponents` member function
    using UnderlyingFunctionType::getNumberOfComponents;
    //! \brief inherit the `getSpace` member function
    using UnderlyingFunctionType::getSpace;
    //! \brief inherit the call operators of the view
    using UnderlyingTensorView::operator();
    //! \brief inherit the `data` member function
    using UnderlyingFunctionType::data;
    //     /*!
    //      * \brief fill the structure using raw data
    //      *
    //      * \param[in] eh: error handler
    //      * \param[in] values: raw data
    //      */
    //     constexpr bool fill(AbstractErrorHandler&, std::span<const real>)
    //     noexcept;
    //     /*!
    //      * \brief fill the structure using raw data
    //      *
    //      * \param[in] eh: error handler
    //      * \param[in] values: raw data
    //      */
    //     constexpr bool fill(AbstractErrorHandler&,
    //                         std::initializer_list<real>) noexcept;
    //! \brief destructor
    constexpr ~TensorialFunction();

   protected:
    /*!
     * \brief make the `check` member function of the view protected
     *
     * Making it inaccessible prevents the `TensorialFunction` class from being
     * treated as an evaluator: a view must be created to evaluate a function.
     */
    using TensorView<Function<Space,
                              tfel::math::getUnderlyingArrayMinimalSize<
                                  typename TensorType::indexing_policy>()>,
                     TensorType>::check;
  };

  /*!
   * \brief return the underlying space
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function
   * \param[in] f: tensorial function
   * \return the space on which the function is defined
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr decltype(auto) getSpace(
      const TensorialFunction<Space, TensorType>&);

  /*!
   * \brief a function whose values are symmetric tensors of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the symmetric tensors, which must be equal to 1, 2
   * or 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using StensorFunction =
      TensorialFunction<Space, tfel::math::stensor<N, real>>;
  /*!
   * \brief a function whose values are tensors of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the tensors, which must be equal to 1, 2 or 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using TensorFunction = TensorialFunction<Space, tfel::math::tensor<N, real>>;
  /*!
   * \brief a function whose values are fourth order tensors mapping symmetric
   * tensors to symmetric tensors, of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the symmetric tensors mapped, which must be equal
   * to 1, 2 or 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using ST2toST2Function =
      TensorialFunction<Space, tfel::math::st2tost2<N, real>>;
  /*!
   * \brief a function whose values are fourth order tensors mapping tensors to
   * symmetric tensors, of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the tensors mapped, which must be equal to 1, 2 or
   * 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using T2toST2Function =
      TensorialFunction<Space, tfel::math::t2tost2<N, real>>;
  /*!
   * \brief a function whose values are fourth order tensors mapping symmetric
   * tensors to tensors, of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the tensors mapped, which must be equal to 1, 2 or
   * 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using ST2toT2Function =
      TensorialFunction<Space, tfel::math::st2tot2<N, real>>;
  /*!
   * \brief a function whose values are fourth order tensors mapping tensors to
   * tensors, of dimension `N`
   *
   * \tparam Space: type of the functional space
   * \tparam N: dimension of the tensors mapped, which must be equal to 1, 2 or
   * 3
   */
  template <FunctionalSpaceConcept Space, unsigned short N>
  using T2toT2Function = TensorialFunction<Space, tfel::math::t2tot2<N, real>>;

  /*!
   * \brief convert a tensorial function to a mutable view
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function
   * \param[in] f: tensorial function
   * \return a mutable view of the values of the function
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(TensorialFunction<Space, TensorType>&);

  /*!
   * \brief convert a tensorial function to a read-only view
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function
   * \param[in] f: tensorial function
   * \return a read-only view of the values of the function
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(const TensorialFunction<Space, TensorType>&);

  /*!
   * \brief creating a view of a temporary tensorial function is forbidden
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  auto view(TensorialFunction<Space, TensorType>&&) = delete;

  /*!
   * \brief return the number of components
   * \tparam Space: type of the functional space
   * \tparam TensorType: type of the tensorial objects returned by the function
   * \param[in] f: tensorial function
   * \return the number of components of the function, which is the number of
   * components of the tensorial type
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorialFunction<Space, TensorType>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorialFunction.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_HXX */
