/*!
 * \file   MGIS/Function/TFEL/Quantity.hxx
 * \brief  File which defines the concepts and the modifiers used to handle the
 * values of a function or an evaluator as TFEL's quantities, i.e. scalar
 * values associated with a physical unit.
 * \author Thomas Helfer
 * \date   16/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_QUANTITY_HXX
#define LIB_MGIS_FUNCTION_TFEL_QUANTITY_HXX

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#include "TFEL/Math/qt.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/QuantityView.hxx"
#include "MGIS/Function/TFEL/QuantityModifier.hxx"

namespace mgis::function::internals {

  /*!
   * \brief specialization for a quantity owning its value, which has one
   * component
   *
   * \tparam UnitType: unit of the quantity
   */
  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::qt<UnitType, real>> {
    //! \brief number of components
    static constexpr size_type value = 1;
  };

  /*!
   * \brief specialization for a mutable reference to a quantity, which has one
   * component
   *
   * \tparam UnitType: unit of the quantity
   */
  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::qt_ref<UnitType, real>> {
    //! \brief number of components
    static constexpr size_type value = 1;
  };

  /*!
   * \brief specialization for a read-only reference to a quantity, which has
   * one component
   *
   * \tparam UnitType: unit of the quantity
   */
  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::const_qt_ref<UnitType, real>> {
    //! \brief number of components
    static constexpr size_type value = 1;
  };

  /*!
   * \brief structure used to modify a function or an evaluator so that it
   * returns its values as quantities, i.e. scalar values associated with a
   * physical unit
   *
   * \tparam UnitType: unit of the quantities returned by the modified function
   * or evaluator
   */
  template <::tfel::math::unit::UnitConcept UnitType>
  struct quantity_modifier {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;
    /*!
     * \brief create a new view
     *
     * \tparam QualifiedFunctionArgumentType: type of the function, which may
     * be a qualified type (reference, const, etc...)
     * \param[in] f: function
     * \return a view returning the values of the function as quantities with
     * the unit `UnitType`
     * \note only available if the number of components of the function is
     * compatible with a quantity, i.e. is equal to one or is only known at run
     * time
     */
    template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType>
    constexpr auto operator()(QualifiedFunctionArgumentType&& f) const
        requires(internals::checkNumberOfComponentsCompatibility<
                 QualifiedFunctionArgumentType,
                 1u>()) {
      return QuantityView<std::decay_t<QualifiedFunctionArgumentType>,
                          UnitType>(
          std::forward<QualifiedFunctionArgumentType>(f));
    }
    /*!
     * \brief create a new modifier
     *
     * \tparam EvaluatorType: evaluator type
     * \param[in] e: evaluator
     * \return an evaluator returning the values computed by the given
     * evaluator as quantities with the unit `UnitType`
     * \note only available if the number of components of the evaluator is
     * compatible with a quantity, i.e. is equal to one or is only known at run
     * time
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType& e) const requires(
        checkEvaluatorNumberOfComponentsCompatibility<EvaluatorType, 1u>()) {
      return QuantityModifier<EvaluatorType, UnitType>(e);
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief concept defining the types which are quantities of real numbers
   *
   * A type matches this concept if, once the qualifiers removed, it is a TFEL
   * quantity, i.e. a value associated with a physical unit, whose underlying
   * numeric type is `real`, whatever the way the value is owned (by value, by
   * mutable reference or by read-only reference)
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept QuantityConcept = (::tfel::math::QuantityConcept<std::decay_t<T>>)&&(
      std::same_as<::tfel::math::base_type<std::decay_t<T>>, real>);

  /*!
   * \brief create a view returning the values of a function as quantities
   *
   * \tparam QualifiedFunctionArgumentType: type of the function, which may be
   * a qualified type (reference, const, etc...)
   * \tparam UnitType: unit of the quantities returned by the view
   * \param[in] f: function
   * \param[in] m: modifier
   * \return a view returning the values of the function as quantities with the
   * unit `UnitType`
   * \note only available if the number of components of the function is
   * compatible with a quantity, i.e. is equal to one or is only known at run
   * time
   */
  template <ViewableFunctionArgumentConcept QualifiedFunctionArgumentType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto operator|(QualifiedFunctionArgumentType&& f,
                           const internals::quantity_modifier<UnitType>& m)  //
      requires(internals::checkNumberOfComponentsCompatibility<
               QualifiedFunctionArgumentType,
               1u>()) {
    return m(std::forward<QualifiedFunctionArgumentType>(f));
  }  // end of operator|

  /*!
   * \brief modifier used to return the values of a function or an evaluator as
   * quantities
   *
   * \tparam UnitType: unit of the quantities
   */
  template <::tfel::math::unit::UnitConcept UnitType>
  inline constexpr auto as_quantity = internals::quantity_modifier<UnitType>{};

  /*!
   * \brief concept defining the types which can address a scalar value, either
   * a real scalar or a quantity of real numbers
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept ScalarConcept = (QuantityConcept<T>) ||
                          (std::same_as<std::decay_t<T>, real>);

  /*!
   * \brief concept defining the non const qualified scalar types, i.e. the
   * scalar types which can be used to modify values
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept MutableScalarConcept = (ScalarConcept<T>) && (!std::is_const_v<T>);

  namespace internals {

    /*!
     * \brief metafunction giving the modifier used to view the values of a
     * function or an evaluator as scalars of type `T`
     *
     * \tparam T: scalar type, which must satisfy the `MutableScalarConcept`
     * concept
     */
    template <MutableScalarConcept T>
    struct ScalarModifier {
      //! \brief modifier returning quantities with the unit of `T`
      using type = quantity_modifier<::tfel::math::quantity_unit<T>>;
    };

    /*!
     * \brief specialization for real scalars, which have no physical unit
     */
    template <>
    struct ScalarModifier<real> {
      //! \brief modifier returning the first component of the values
      using type = fixed_size_modifier<1>;
    };

  }  // end of namespace internals

  /*!
   * \brief modifier used to return the values of a function or an evaluator as
   * scalars of type `T`
   *
   * If `T` is a quantity, this modifier is `as_quantity` with the unit of `T`.
   * If `T` is `real`, which is the default, this modifier is `as_scalar`
   *
   * \tparam T: scalar type, which must satisfy the `MutableScalarConcept`
   * concept
   */
  template <MutableScalarConcept T = real>
  inline constexpr auto as_qt = typename internals::ScalarModifier<T>::type{};

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITY_HXX */
