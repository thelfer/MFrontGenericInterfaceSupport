/*!
 * \file   MGIS/Function/TFEL/Tensors.hxx
 * \brief  File which defines the modifiers turning scalar valued functions and
 * evaluators into tensor valued ones (`as_stensor`, `as_tensor`, ...) together
 * with the tensorial operations which can be applied to evaluators (`trace`,
 * `det`, `eigen_values`, `rotate`, `logarithmic_strain`, ...) and the concepts
 * characterising evaluators returning the various kinds of tensorial objects.
 * \author Thomas Helfer
 * \date   11/05/2025
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

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORS_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORS_HXX

#include "TFEL/Math/Forward/General.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"
#include "MGIS/Function/TFEL/TensorView.hxx"
#include "MGIS/Function/TFEL/TensorModifier.hxx"
#include "MGIS/Function/TFEL/CoalescedMemoryAccessTensorView.hxx"
#include "MGIS/Function/TFEL/CoalescedMemoryAccessCompositeTensorsView.hxx"
#include "MGIS/Function/TFEL/StridedCoalescedMemoryAccessTensorView.hxx"
#include "MGIS/Function/TFEL/StridedCoalescedMemoryAccessCompositeTensorsView.hxx"

namespace mgis::function::internals {

  /*!
   * \brief partial specialization stating that a mutable reference to a
   * tensorial object is a valid result type for a function
   *
   * \tparam TensorType: type of the tensorial object, which must match the
   * `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct FunctionResultTypeTraits<TensorType&> {
    //! \brief the traits are specialized for this type
    static constexpr auto is_specialized = true;
  };

  /*!
   * \brief partial specialization stating that a view on a tensorial object is
   * a valid result type for a function
   *
   * \tparam TensorType: type of the tensorial object viewed, which must match
   * the `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct FunctionResultTypeTraits<tfel::math::View<TensorType>> {
    //! \brief the traits are specialized for this type
    static constexpr auto is_specialized = true;
  };

  /*!
   * \brief modifier exposing the values computed by an evaluator, or the values
   * of a function, as a tensorial object of type `TensorType`
   *
   * This is the modifier behind the `as_stensor`, `as_tensor`, ... variables
   * defined below. It is stateless, which allows those variables to be defined
   * as `constexpr` global variables.
   *
   * \tparam TensorType: type of the tensorial object to be exposed, which must
   * match the `TensorConcept` concept
   */
  template <TensorConcept TensorType>
  struct tensor_modifier {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;

    /*!
     * \brief expose the values of a function argument as a tensorial object
     * \tparam QualifiedFunctionArgumentType: type of the function argument,
     * which may be a qualified type
     * \param[in] f: function argument whose values are exposed
     * \return a `TensorView` of type `TensorType` viewing the values of the
     * function argument
     * \note only available if the number of components of the function argument
     * is compatible with that of `TensorType` (the check always succeeds when
     * the number of components is only known at run time)
     */
    template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType>
    constexpr auto operator()(QualifiedFunctionArgumentType&&) const
        requires(internals::checkNumberOfComponentsCompatibility<
                 QualifiedFunctionArgumentType,
                 compile_time_size<TensorType>>());
    /*!
     * \brief expose the values computed by an evaluator as a tensorial object
     * \tparam EvaluatorType: type of the underlying evaluator
     * \param[in] e: evaluator whose values are exposed
     * \return a `TensorModifier` returning the values computed by the evaluator
     * as objects of type `TensorType`
     * \note only available if the requirements of the `TensorModifier` class
     * are satisfied for `TensorType` and the given evaluator type
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType&) const requires(
        (areTensorModifierRequirementsSatisfied<TensorType,
                                                std::decay_t<EvaluatorType>>));
  };

  /*!
   * \brief modifier changing the basis of the tensorial values computed by an
   * evaluator, in the direction for which the given rotation matrix is defined
   *
   * The rotation matrix can be given either as a `tfel::math::tmatrix`, as
   * another evaluator providing the rotation matrix, or directly together with
   * the evaluator to rotate, through the associated `operator|`.
   */
  struct RotateModifier {
    //! \brief type of the callable object performing the rotation
    using CallableType = decltype(
        []<typename TensorType>(const TensorType& t,
                                const tfel::math::tmatrix<3, 3, real>& R)  //
        requires((tfel::math::TensorConcept<TensorType>) ||
                 (tfel::math::StensorConcept<TensorType>) ||
                 (tfel::math::ST2toST2Concept<TensorType>)) {
          return tfel::math::change_basis(t, R);
        });

    /*!
     * \brief create a modifier rotating the values computed by an evaluator
     * into the basis given by a rotation matrix
     * \param[in] R: rotation matrix, whose columns are the basis vectors
     * expressed in the initial basis
     * \return a unary operation modifier computing
     * `tfel::math::change_basis(t, R)`
     */
    constexpr auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [R]<typename TensorType>(const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>) ||
                   (tfel::math::ST2toST2Concept<TensorType>)) {
        return tfel::math::change_basis(t, R);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    /*!
     * \brief partially apply the rotation to the evaluator of the second
     * argument, i.e. the evaluator providing the rotation matrix
     * \tparam SecondEvaluatorType: type of the evaluator providing the rotation
     * matrix, which may be a qualified type
     * \param[in] e2: evaluator used to compute the rotation matrix
     * \return a `BinaryOperatorCurrying2` object, which can be applied to the
     * evaluator of the first argument using `operator|`
     */
    template <typename SecondEvaluatorType>
    constexpr auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    /*!
     * \brief build the evaluator of the rotated tensorial values
     * \tparam FirstEvaluatorType: type of the evaluator providing the tensorial
     * values to rotate, which may be a qualified type
     * \tparam SecondEvaluatorType: type of the evaluator providing the rotation
     * matrix, which may be a qualified type
     * \param[in] e1: evaluator used to compute the tensorial values to rotate
     * \param[in] e2: evaluator used to compute the rotation matrix
     * \return a `BinaryOperationModifier2` object returning the rotated values
     */
    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    constexpr auto operator()(FirstEvaluatorType&& e1,
                              SecondEvaluatorType&& e2) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                 (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                 (std::invocable<
                     CallableType,
                     evaluator_result<std::decay_t<FirstEvaluatorType>>,
                     evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
      return BinaryOperationModifier2<CallableType,
                                      std::decay_t<FirstEvaluatorType>,
                                      std::decay_t<SecondEvaluatorType>>(
          std::forward<FirstEvaluatorType>(e1),
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()
  };

  /*!
   * \brief modifier changing the basis of the tensorial values computed by an
   * evaluator in the backwards direction, i.e. using the transpose of the given
   * rotation matrix
   *
   * This is the counterpart of `RotateModifier`: where the latter gives the
   * components of the values in the basis whose vectors are the columns of the
   * rotation matrix, the former gives their components in the initial basis
   * knowing their components in that basis.
   */
  struct RotateBackwardsModifier {
    //! \brief type of the callable object performing the rotation
    using CallableType = decltype(
        []<typename TensorType>(const TensorType& t,
                                const tfel::math::tmatrix<3, 3, real>& R)  //
        requires((tfel::math::TensorConcept<TensorType>) ||
                 (tfel::math::StensorConcept<TensorType>) ||
                 (tfel::math::ST2toST2Concept<TensorType>)) {
          return tfel::math::change_basis(t, tfel::math::transpose(R));
        });

    /*!
     * \brief create a modifier rotating the values computed by an evaluator
     * using the transpose of a rotation matrix
     * \param[in] R: rotation matrix, whose columns are the basis vectors
     * expressed in the initial basis; its transpose is used
     * \return a unary operation modifier computing
     * `tfel::math::change_basis(t, tfel::math::transpose(R))`
     */
    constexpr auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [Rb = tfel::math::transpose(R)]<typename TensorType>(
          const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>) ||
                   (tfel::math::ST2toST2Concept<TensorType>)) {
        return tfel::math::change_basis(t, Rb);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    /*!
     * \brief partially apply the rotation to the evaluator of the second
     * argument, i.e. the evaluator providing the rotation matrix
     * \tparam SecondEvaluatorType: type of the evaluator providing the rotation
     * matrix, which may be a qualified type
     * \param[in] e2: evaluator used to compute the rotation matrix
     * \return a `BinaryOperatorCurrying2` object, which can be applied to the
     * evaluator of the first argument using `operator|`
     */
    template <typename SecondEvaluatorType>
    constexpr auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    /*!
     * \brief build the evaluator of the rotated tensorial values
     * \tparam FirstEvaluatorType: type of the evaluator providing the tensorial
     * values to rotate, which may be a qualified type
     * \tparam SecondEvaluatorType: type of the evaluator providing the rotation
     * matrix, which may be a qualified type
     * \param[in] e1: evaluator used to compute the tensorial values to rotate
     * \param[in] e2: evaluator used to compute the rotation matrix
     * \return a `BinaryOperationModifier2` object returning the rotated values
     */
    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    constexpr auto operator()(FirstEvaluatorType&& e1,
                              SecondEvaluatorType&& e2) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                 (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                 (std::invocable<
                     CallableType,
                     evaluator_result<std::decay_t<FirstEvaluatorType>>,
                     evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
      return BinaryOperationModifier2<CallableType,
                                      std::decay_t<FirstEvaluatorType>,
                                      std::decay_t<SecondEvaluatorType>>(
          std::forward<FirstEvaluatorType>(e1),
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()
  };   // end of RotateBackwardModifier

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief concept satisfied by evaluators returning a symmetric tensor of
   * dimension `1` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Stensor1DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::StensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a symmetric tensor of
   * dimension `2` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Stensor2DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::StensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a symmetric tensor of
   * dimension `3` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Stensor3DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::StensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a symmetric tensor of any
   * dimension whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept StensorEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::StensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning an unsymmetric tensor of
   * dimension `1` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Tensor1DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::TensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning an unsymmetric tensor of
   * dimension `2` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Tensor2DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::TensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning an unsymmetric tensor of
   * dimension `3` whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept Tensor3DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::TensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning an unsymmetric tensor of
   * any dimension whose base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept TensorEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::TensorConcept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to symmetric tensors, of dimension `1`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toST21DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                 //
      (::tfel::math::ST2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to symmetric tensors, of dimension `2`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toST22DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                 //
      (::tfel::math::ST2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to symmetric tensors, of dimension `3`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toST23DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                 //
      (::tfel::math::ST2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to symmetric tensors, of any dimension, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toST2EvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::ST2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to unsymmetric tensors, of dimension `1`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toT21DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::ST2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to unsymmetric tensors, of dimension `2`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toT22DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::ST2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to unsymmetric tensors, of dimension `3`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toT23DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::ST2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping symmetric tensors to unsymmetric tensors, of any dimension, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept ST2toT2EvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::ST2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to symmetric tensors, of dimension `1`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toST21DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::T2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to symmetric tensors, of dimension `2`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toST22DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::T2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to symmetric tensors, of dimension `3`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toST23DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                                //
      (::tfel::math::T2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to symmetric tensors, of any dimension, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toST2EvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::T2toST2Concept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to unsymmetric tensors, of dimension `1`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toT21DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::T2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       1u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to unsymmetric tensors, of dimension `2`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toT22DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::T2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       2u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to unsymmetric tensors, of dimension `3`, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toT23DEvaluatorConcept =
      (EvaluatorConcept<EvaluatorType>)&&                               //
      (::tfel::math::T2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (::tfel::math::space_dimension<evaluator_result<EvaluatorType>> ==
       3u) &&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor
   * mapping unsymmetric tensors to unsymmetric tensors, of any dimension, whose
   * base numeric type is `real`
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept T2toT2EvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&  //
      (::tfel::math::T2toT2Concept<evaluator_result<EvaluatorType>>)&&  //
      (std::same_as<
          std::decay_t<::tfel::math::base_type<
              ::tfel::math::numeric_type<evaluator_result<EvaluatorType>>>>,
          real>);

  /*!
   * \brief concept satisfied by evaluators returning a fourth order tensor of
   * any kind, i.e. matching one of the `ST2toST2EvaluatorConcept`,
   * `T2toST2EvaluatorConcept`, `ST2toT2EvaluatorConcept` or
   * `T2toT2EvaluatorConcept` concepts
   *
   * \tparam EvaluatorType: type of the evaluator considered
   */
  template <typename EvaluatorType>
  concept FourthOrderTensorEvaluatorConcept =
      ST2toST2EvaluatorConcept<EvaluatorType> ||
      T2toST2EvaluatorConcept<EvaluatorType> ||
      ST2toT2EvaluatorConcept<EvaluatorType> ||
      T2toT2EvaluatorConcept<EvaluatorType>;

  /*!
   * \brief expose the values of a function argument as a tensorial object
   *
   * This is the entry point used by expressions such as `f | as_stensor<3>`.
   *
   * \tparam QualifiedFunctionArgumentType: type of the function argument, which
   * may be a qualified type
   * \tparam TensorType: type of the tensorial object to be exposed
   * \param[in] f: function argument whose values are exposed
   * \param[in] m: modifier stating the type of the tensorial object to be
   * exposed
   * \return a `TensorView` viewing the values of the function argument as
   * objects of type `TensorType`
   * \note only available if the number of components of the function argument
   * is compatible with that of `TensorType`
   */
  template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType,
            TensorConcept TensorType>
  constexpr auto operator|(QualifiedFunctionArgumentType&&,
                           const internals::tensor_modifier<TensorType>&)  //
      requires(internals::checkNumberOfComponentsCompatibility<
               QualifiedFunctionArgumentType,
               compile_time_size<TensorType>>());

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * symmetric tensors of dimension `N`
   *
   * \tparam N: dimension of the symmetric tensors, which must be 1, 2 or 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_stensor =
          internals::tensor_modifier<tfel::math::stensor<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * unsymmetric tensors of dimension `N`
   *
   * \tparam N: dimension of the tensors, which must be 1, 2 or 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_tensor =
          internals::tensor_modifier<tfel::math::tensor<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * finite size arrays of `N` values
   *
   * \tparam N: number of values stored in the array
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  inline constexpr auto as_fsarray =
      internals::tensor_modifier<tfel::math::fsarray<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * tiny vectors of `N` components
   *
   * \tparam N: number of components of the vector
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  inline constexpr auto as_tvector =
      internals::tensor_modifier<tfel::math::tvector<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * tiny matrices with `N` rows and `M` columns
   *
   * \tparam N: number of rows of the matrix
   * \tparam M: number of columns of the matrix
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N,
            unsigned short M,
            TensorValueConcept ValueType = real>
  inline constexpr auto as_tmatrix =
      internals::tensor_modifier<tfel::math::tmatrix<N, M, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * fourth order tensors mapping symmetric tensors to symmetric tensors, of
   * dimension `N`
   *
   * \tparam N: dimension of the symmetric tensors mapped, which must be 1, 2 or
   * 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_st2tost2 =
          internals::tensor_modifier<tfel::math::st2tost2<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * fourth order tensors mapping unsymmetric tensors to symmetric tensors, of
   * dimension `N`
   *
   * \tparam N: dimension of the tensors mapped, which must be 1, 2 or 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_t2tost2 =
          internals::tensor_modifier<tfel::math::t2tost2<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * fourth order tensors mapping symmetric tensors to unsymmetric tensors, of
   * dimension `N`
   *
   * \tparam N: dimension of the tensors mapped, which must be 1, 2 or 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_st2tot2 =
          internals::tensor_modifier<tfel::math::st2tot2<N, ValueType>>{};

  /*!
   * \brief modifier exposing the values of a function or of an evaluator as
   * fourth order tensors mapping unsymmetric tensors to unsymmetric tensors, of
   * dimension `N`
   *
   * \tparam N: dimension of the tensors mapped, which must be 1, 2 or 3
   * \tparam ValueType: type of the values, which must match the
   * `TensorValueConcept` concept and defaults to `real`
   */
  template <unsigned short N, TensorValueConcept ValueType = real>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_t2tot2 =
          internals::tensor_modifier<tfel::math::t2tot2<N, ValueType>>{};

  /*!
   * \brief modifier computing the trace of the tensorial values computed by an
   * evaluator
   *
   * The values must be symmetric or unsymmetric tensors; the modifier returns a
   * scalar valued result.
   */
  inline constexpr auto trace = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::trace(t);
      });

  /*!
   * \brief modifier computing the determinant of the tensorial values computed
   * by an evaluator
   *
   * The values must be symmetric or unsymmetric tensors; the modifier returns a
   * scalar valued result.
   */
  inline constexpr auto det = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::det(t);
      });

  /*!
   * \brief modifier computing the eigen values of the symmetric tensorial
   * values computed by an evaluator
   *
   * The values must be symmetric tensors; the modifier returns a vector holding
   * as many eigen values as the dimension of the tensors.
   *
   * \tparam esolver: eigen solver used to compute the eigen values, which
   * defaults to `tfel::math::stensor_common::TFELEIGENSOLVER`
   */
  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto eigen_values = internals::unary_operation_modifier2(
      [](const tfel::math::StensorConcept auto& s) {
        auto s2 = tfel::math::EvaluationResult<decltype(s)>(s);
        return s2.template computeEigenValues<esolver>();
      });

  //! \brief modifier changing the basis of the tensorial values of an evaluator
  inline constexpr auto rotate = internals::RotateModifier{};

  //! \brief counterpart of `rotate`, using the transpose of the rotation matrix
  inline constexpr auto rotate_backwards = internals::RotateBackwardsModifier{};

  /*!
   * \brief modifier computing the logarithmic strain from the deformation
   * gradients computed by an evaluator
   *
   * The values must be unsymmetric tensors, interpreted as deformation
   * gradients `F`. The logarithmic strain is obtained from the Green-Lagrange
   * strain tensor `e = (F^T F - I) / 2`: its eigen values `vp` are transformed
   * into `log(1 + 2 vp) / 2` and the result is rebuilt from the transformed
   * eigen values and the eigen vectors of `e`.
   *
   * \tparam esolver: eigen solver used to compute the eigen vectors, which
   * defaults to `tfel::math::stensor_common::TFELEIGENSOLVER`
   */
  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto logarithmic_strain =
      internals::unary_operation_modifier2(
          [](const tfel::math::TensorConcept auto& F) {
            using namespace tfel::math;
            constexpr auto N = getSpaceDimension<decltype(F)>();
            const auto e = computeGreenLagrangeTensor(F);
            const auto [vp, m] = e.template computeEigenVectors<esolver>();
            const auto log_vp =
                map([](const real x) { return std::log1p(2 * x) / 2; }, vp);
            return stensor<N, real>::computeIsotropicFunction(log_vp, m);
          });

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/Tensors.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORS_HXX */
