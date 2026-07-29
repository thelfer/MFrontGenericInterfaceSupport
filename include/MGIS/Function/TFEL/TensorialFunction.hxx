/*!
 * \file   MGIS/Function/TFEL/TensorialFunction.hxx
 * \brief
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
#include "MGIS/Function/TFEL/Tensors.hxx"

namespace mgis::function {

  /*!
   * \brief default implementation of a function
   *
   * \tparam Space: type of the functional space
   * \tparam TensorType: tensor type
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
    //! \brief a simple alias
    using UnderlyingFunctionType =
        Function<Space,
                 tfel::math::getUnderlyingArrayMinimalSize<
                     typename TensorType::indexing_policy>()>;
    //! \brief a simple alias
    using UnderlyingTensorView =
        TensorView<Function<Space,
                            tfel::math::getUnderlyingArrayMinimalSize<
                                typename TensorType::indexing_policy>()>,
                   TensorType>;
    /*!
     * \brief constructor from a space and a data size
     * \param[in] eh: error handler
     * \param[in] s: space
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&, const Space&);
    /*!
     * \brief constructor from a space
     * Vieparam[in] s: space
     */
    constexpr TensorialFunction(const Space&);
    /*!
     * \brief constructor from a space
     * \param[in] s: space
     */
    template <bool doPreconditionsCheck>
    constexpr TensorialFunction(const PreconditionsCheck<doPreconditionsCheck>&,
                                const Space&);
    //! \brief move constructor
    constexpr TensorialFunction(TensorialFunction&&);
    //! \brief copy constructor
    constexpr TensorialFunction(const TensorialFunction&);
    //! \brief return a view of the function
    constexpr UnderlyingTensorView& view() noexcept;
    //! \brief return a view of the function
    constexpr const UnderlyingTensorView& view() const noexcept;
    //
    using UnderlyingFunctionType::getNumberOfComponents;
    using UnderlyingFunctionType::getSpace;
    using UnderlyingTensorView::operator();
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
    /*
     * This function is made protected to avoid Function from being treated
     * as an evaluator
     */
    using TensorView<Function<Space,
                              tfel::math::getUnderlyingArrayMinimalSize<
                                  typename TensorType::indexing_policy>()>,
                     TensorType>::check;
  };

  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr decltype(auto) getSpace(
      const TensorialFunction<Space, TensorType>&);

  template <FunctionalSpaceConcept Space, unsigned short N>
  using StensorFunction =
      TensorialFunction<Space, tfel::math::stensor<N, real>>;
  template <FunctionalSpaceConcept Space, unsigned short N>
  using TensorFunction = TensorialFunction<Space, tfel::math::tensor<N, real>>;
  template <FunctionalSpaceConcept Space, unsigned short N>
  using ST2toST2Function =
      TensorialFunction<Space, tfel::math::st2tost2<N, real>>;
  template <FunctionalSpaceConcept Space, unsigned short N>
  using T2toST2Function =
      TensorialFunction<Space, tfel::math::t2tost2<N, real>>;
  template <FunctionalSpaceConcept Space, unsigned short N>
  using ST2toT2Function =
      TensorialFunction<Space, tfel::math::st2tot2<N, real>>;
  template <FunctionalSpaceConcept Space, unsigned short N>
  using T2toT2Function = TensorialFunction<Space, tfel::math::t2tot2<N, real>>;

  /*!
   * \brief convert a tensorial function to a immutable view
   * \param[in] f: function
   */
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr auto view(const TensorialFunction<Space, TensorType>&);

  //! \return the number of components
  template <FunctionalSpaceConcept Space, TensorConcept TensorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorialFunction<Space, TensorType>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorialFunction.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORIALFUNCTION_HXX */
