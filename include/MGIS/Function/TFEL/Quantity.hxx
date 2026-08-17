/*!
 * \file   MGIS/Function/TFEL/Quantity.hxx
 * \brief
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
#include "MGIS/Function/TFEL/QuantityView.hxx"
#include "MGIS/Function/TFEL/QuantityModifier.hxx"

namespace mgis::function::internals {

  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::qt<UnitType, real>> {
    static constexpr size_type value = 1;
  };

  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::qt_ref<UnitType, real>> {
    static constexpr size_type value = 1;
  };

  template <::tfel::math::unit::UnitConcept UnitType>
  struct CompileTimeSize<::tfel::math::const_qt_ref<UnitType, real>> {
    static constexpr size_type value = 1;
  };

  template <::tfel::math::unit::UnitConcept UnitType>
  struct quantity_modifier {
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;
    /*!
     * \brief create a new view
     * \param[in] f: function type
     */
    template <FunctionConcept FunctionType>
    constexpr auto operator()(FunctionType& f) const
        requires(number_of_components<FunctionType> == dynamic_extent
                     ? true
                     : 1 == number_of_components<FunctionType>) {
      return QuantityView<FunctionType, UnitType>(f);
    }
    /*!
     * \brief create a new modifier
     * \param[in] e: evaluator type
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType& e) const
        requires(number_of_components<EvaluatorType> == dynamic_extent
                     ? true
                     : 1 == number_of_components<EvaluatorType>) {
      return QuantityModifier<EvaluatorType, UnitType>(e);
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <typename T>
  concept QuantityConcept = (::tfel::math::QuantityConcept<std::decay_t<T>>)&&(
      std::same_as<::tfel::math::base_type<std::decay_t<T>>, real>);

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto operator|(FunctionType& f,
                           const internals::quantity_modifier<UnitType>& m)  //
      requires(number_of_components<FunctionType> == dynamic_extent
                   ? true
                   : 1 == number_of_components<FunctionType>) {
    return m(f);
  }  // end of operator|

  template <::tfel::math::unit::UnitConcept UnitType>
  inline constexpr auto as_quantity = internals::quantity_modifier<UnitType>{};

  template <typename T>
  concept ScalarConcept = (QuantityConcept<T>) ||
                          (std::same_as<std::decay_t<T>, real>);

  template <typename T>
  concept MutableScalarConcept = (ScalarConcept<T>) || (!std::is_const_v<T>);

  namespace internals {

    template <MutableScalarConcept T>
    struct ScalarModifier{
      using type = quantity_modifier<::tfel::math::quantity_unit<T>>;
    };

    template <>
    struct ScalarModifier<real> {
      using type = fixed_size_modifier<1>;
    };

  }  // end of namespace internals

  template <MutableScalarConcept T = real>
  inline constexpr auto as_qt =
      typename internals::ScalarModifier<T>::type{};

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITY_HXX */
