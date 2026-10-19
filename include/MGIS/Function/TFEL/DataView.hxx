/*!
 * \file   MGIS/Function/TFEL/DataView.hxx
 * \brief  This file specializes the `IsDataView` class for TFEL' types
 * \author Thomas Helfer
 * \date   02/10/2026
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX

#include "TFEL/Math/Forward/qt.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "TFEL/Math/Array/CoalescedView.hxx"
#include "MGIS/Function/DataViewConcept.hxx"

namespace mgis::function::internals {

  //! \brief partial specialisation for qt_ref
  template <::tfel::math::unit::UnitConcept UnitType,
            ::tfel::math::StandardArithmeticTypeConcept ValueType>
  struct IsDataView<::tfel::math::qt_ref<UnitType, ValueType>>
      : std::true_type {};

  //! \brief partial specialisation for const_qt_ref
  template <::tfel::math::unit::UnitConcept UnitType,
            ::tfel::math::StandardArithmeticTypeConcept ValueType>
  struct IsDataView<::tfel::math::const_qt_ref<UnitType, ValueType>>
      : std::true_type {};

  //! \brief partial specialisation for tfel::math::View
  template <typename MappedType, typename IndexingPolicyType>
  struct IsDataView<::tfel::math::View<MappedType, IndexingPolicyType>>
    : std::true_type{};

  //! \brief partial specialisation for tfel::math::CoalescedView
  template <
      ::tfel::math::MappableMathObjectUsingCoalescedViewConcept MappedType,
      typename IndexingPolicyType>
  struct IsDataView<::tfel::math::CoalescedView<MappedType, IndexingPolicyType>>
    : std::true_type{};

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX */
