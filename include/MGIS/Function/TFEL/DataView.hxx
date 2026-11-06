/*!
 * \file   MGIS/Function/TFEL/DataView.hxx
 * \brief  This file specializes the `IsDataView` class for TFEL's types,
 *         which all view values owned by other objects
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

  /*!
   * \brief partial specialisation for `tfel::math::qt_ref`, a quantity which
   * wraps a mutable reference to a value owned by another object
   *
   * \tparam UnitType: unit of the quantity
   * \tparam ValueType: numeric type of the referred value
   */
  template <::tfel::math::unit::UnitConcept UnitType,
            ::tfel::math::StandardArithmeticTypeConcept ValueType>
  struct IsDataView<::tfel::math::qt_ref<UnitType, ValueType>>
      : std::true_type {};

  /*!
   * \brief partial specialisation for `tfel::math::const_qt_ref`, a quantity
   * which wraps a read-only reference to a value owned by another object
   *
   * \tparam UnitType: unit of the quantity
   * \tparam ValueType: numeric type of the referred value
   */
  template <::tfel::math::unit::UnitConcept UnitType,
            ::tfel::math::StandardArithmeticTypeConcept ValueType>
  struct IsDataView<::tfel::math::const_qt_ref<UnitType, ValueType>>
      : std::true_type {};

  /*!
   * \brief partial specialisation for `tfel::math::View`, a view of a
   * mathematical object which may use a non contiguous indexing policy
   *
   * \tparam MappedType: type of the mathematical object viewed
   * \tparam IndexingPolicyType: policy used to address the values of the
   * viewed object
   */
  template <typename MappedType, typename IndexingPolicyType>
  struct IsDataView<::tfel::math::View<MappedType, IndexingPolicyType>>
      : std::true_type {};

  /*!
   * \brief partial specialisation for `tfel::math::CoalescedView`, a view of
   * a mathematical object whose values are stored in a coalesced memory
   *
   * \tparam MappedType: type of the mathematical object viewed, which must be
   * mappable using a coalesced view
   * \tparam IndexingPolicyType: policy used to address the values of the
   * viewed object
   */
  template <
      ::tfel::math::MappableMathObjectUsingCoalescedViewConcept MappedType,
      typename IndexingPolicyType>
  struct IsDataView<::tfel::math::CoalescedView<MappedType, IndexingPolicyType>>
      : std::true_type {};

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX */
