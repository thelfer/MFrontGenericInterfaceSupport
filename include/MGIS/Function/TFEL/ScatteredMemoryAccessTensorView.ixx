/*!
 * \file   MGIS/Function/TFEL/ScatteredMemoryAccessTensorView.ixx
 * \brief  File which contains the implementation of the
 * `ScatteredMemoryAccessTensorView` class declared in
 * `MGIS/Function/TFEL/ScatteredMemoryAccessTensorView.hxx`.
 * \author Thomas Helfer
 * \date   27/10/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_SCATTEREDMEMORYACCESSTENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_SCATTEREDMEMORYACCESSTENSORVIEW_IXX

namespace mgis::function {

  // the call operators build a scattered view from the array of the addresses
  // of the values of the components at the given integration point, returned by
  // the `getValuesPointers` member function of the base class
  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename ScatteredMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::MutableValues
  ScatteredMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type
          i) requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                      (!hasElementWorkspace<Space>)) {
    return tfel::math::ScatteredView<TensorType>(this->getValuesPointers(i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename ScatteredMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::MutableValues
  ScatteredMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type e,
      const size_type q)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return tfel::math::ScatteredView<TensorType>(this->getValuesPointers(e, q));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename ScatteredMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::ConstValues
  ScatteredMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type i) const  //
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return tfel::math::ScatteredView<const TensorType>(
        this->getValuesPointers(i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr typename ScatteredMemoryAccessTensorView<Space,
                                                     TensorType,
                                                     is_mutable>::ConstValues
  ScatteredMemoryAccessTensorView<Space, TensorType, is_mutable>::operator()(
      const size_type e, const size_type q) const  //
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return tfel::math::ScatteredView<const TensorType>(
        this->getValuesPointers(e, q));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr const Space& getSpace(
      const ScatteredMemoryAccessTensorView<Space, TensorType, is_mutable>&
          v) MGIS_NOEXCEPT {
    return v.getSpace();
  }  // end of getSpace

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr bool check(
      AbstractErrorHandler&,
      const ScatteredMemoryAccessTensorView<Space,
                                            TensorType,
                                            is_mutable>&) MGIS_NOEXCEPT {
    return true;
  }  // end of check

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr size_type getNumberOfComponents(
      const ScatteredMemoryAccessTensorView<Space,
                                            TensorType,
                                            is_mutable>&) MGIS_NOEXCEPT {
    return compile_time_size<TensorType>;
  }  // end of getNumberOfComponents

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_SCATTEREDMEMORYACCESSTENSORVIEW_IXX */
