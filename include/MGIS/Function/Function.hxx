/*!
 * \file   MGIS/Function/Function.hxx
 * \brief File which defines the `FunctionView` class, which describes a
 * function whose values are stored in an external memory region, the
 * `Function` class, which owns the values of a function, and the helper
 * functions used to query them and to create views.
 * \author Thomas Helfer
 * \date   11/06/2020
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTION_HXX
#define LIB_MGIS_FUNCTION_FUNCTION_HXX

#include <span>
#include <limits>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  /*!
   * \brief class describing the size (number of components) of the data
   * manipulated by a function, when this size is known at compile time
   *
   * \tparam data_size: number of components of the function values, which
   * must be greater than zero
   */
  template <size_type data_size>
  requires(data_size > 0) struct FunctionDataSize {
    //! \return true if the function is a scalar valued function
    [[nodiscard]] constexpr bool isScalar() const MGIS_NOEXCEPT;
    //! \return the number of components of the function values
    [[nodiscard]] constexpr size_type getNumberOfComponents() const
        MGIS_NOEXCEPT;
  };

  /*!
   * \brief specialisation of `FunctionDataSize` used when the number of
   * components of the function is only known at run time
   *
   * The number of components is stored in the protected data member
   * `data_size`, which is set by the constructors of the derived classes.
   */
  template <>
  struct FunctionDataSize<dynamic_extent> {
    //! \return true if the function is a scalar valued function
    [[nodiscard]] constexpr bool isScalar() const MGIS_NOEXCEPT;
    //! \return the number of components of the function values
    [[nodiscard]] constexpr size_type getNumberOfComponents() const
        MGIS_NOEXCEPT;

   protected:
    //! \brief number of components of the function values
    size_type data_size = size_type{};
  };  // end of FunctionDataSize<dynamic_extent>

  /*!
   * \brief class describing the stride of the data manipulated by a
   * function, i.e. the distance between the values associated with two
   * successive integration points, when this stride is known at compile time
   *
   * \tparam data_stride: stride of the function data, which must be
   * greater than zero
   */
  template <size_type data_stride>
  requires(data_stride > 0) struct FunctionDataStride {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    [[nodiscard]] constexpr size_type getDataStride() const MGIS_NOEXCEPT;
  };

  /*!
   * \brief specialisation of `FunctionDataStride` used when the stride is
   * only known at run time
   *
   * The stride is stored in the protected data member `data_stride`, which
   * is set by the constructors of the derived classes.
   */
  template <>
  struct FunctionDataStride<dynamic_extent> {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    [[nodiscard]] constexpr size_type getDataStride() const MGIS_NOEXCEPT;

   protected:
    //! \brief stride of the function data
    size_type data_stride = size_type{};
  };  // end of FunctionDataStride<dynamic_extent>

  /*!
   * \brief a simple data structure describing how the data of a function is
   * mapped in memory, which is used as a non-type template argument
   *
   * The default value, whose data size and data stride are both equal to
   * `dynamic_extent`, corresponds to a layout which is entirely determined
   * at run time.
   */
  struct FunctionDataLayoutDescription {
    //! \brief number of components, `dynamic_extent` if unknown at run time
    size_type data_size = dynamic_extent;
    //! \brief stride, `dynamic_extent` if only known at run time
    size_type data_stride = dynamic_extent;
  };

  /*!
   * \brief checks if the given data layout description has properties which
   * are only known at run time
   * \param[in] layout: data layout description
   * \return true if the data size or the data stride is `dynamic_extent`
   */
  [[nodiscard]] constexpr bool has_dynamic_properties(
      const FunctionDataLayoutDescription&);

  /*!
   * \brief a simple data structure describing how the data of a
   * function is mapped in memory
   *
   * This class gathers the number of components and the stride of the data
   * of a function through the `FunctionDataSize` and `FunctionDataStride`
   * base classes, which store those properties at compile time or at run
   * time depending on the value of `layout`.
   *
   * \tparam layout: data layout description. The `requires` clause states
   * that the data size and the data stride are greater than zero and that
   * the data size is not greater than the data stride when both are known
   * at compile time.
   */
  template <FunctionDataLayoutDescription layout>
    requires((layout.data_size > 0) && (layout.data_stride > 0) &&
             ((layout.data_size == dynamic_extent) ||
              (layout.data_stride == dynamic_extent) ||
              (layout.data_size <= layout.data_stride)))
  struct FunctionDataLayout : public FunctionDataSize<layout.data_size>,
                              public FunctionDataStride<layout.data_stride> {
    //! \brief default constructor
    FunctionDataLayout() = default;
    //! \brief move constructor
    FunctionDataLayout(FunctionDataLayout&&) = default;
    //! \brief copy constructor
    FunctionDataLayout(const FunctionDataLayout&) = default;
    //! \brief move assignment operator
    FunctionDataLayout& operator=(FunctionDataLayout&&) = default;
    //! \brief copy assignment operator
    FunctionDataLayout& operator=(const FunctionDataLayout&) = default;
    //! \brief destructor
    ~FunctionDataLayout() = default;

   protected:
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] i: integration point
     */
    [[nodiscard]] constexpr size_type getDataOffset(const size_type) const
        MGIS_NOEXCEPT;
  };  // end of struct FunctionDataLayout

  /*!
   * \brief function defined on a space, whose values are stored in an
   * external memory region.
   *
   * The `FunctionView` defines a function using an external memory region.
   *
   * This memory region may contain more data than the one associated with the
   * function as illustrated by the following figure:
   *
   * |---------------------------------------------------------------|
   * <-                         Raw data                            ->
   * |---------------------------------------|
   * <- Data of the first element          -->
   * |---------------|                       | <-function data->
   *                 ^                       ^
   *                 |                       |
   *             data_size                   |
   *                                     data_stride
   *
   * The size of all the data (including the one not related to the
   * function) associated with one element of the space is called the
   * `data_stride` in the `FunctionView` class.
   *
   * The size of the data hold by the function per element of the space,
   * i.e. the number of components of the function is given by
   * `data_size`.
   *
   * The values of the view can be accessed either with the offset of an
   * integration point, if the space matches the `LinearElementSpaceConcept`
   * concept, or with the index of a cell and the index of a quadrature
   * point, if the space matches the `LinearQuadratureSpaceConcept` concept.
   *
   * A `FunctionView` is a lightweight view: it can be copied at a
   * negligible cost, as stated by the specialisation of
   * `LightweightViewTraits`, and it matches the `EvaluatorConcept` concept:
   * it can be used as an evaluator, which returns constant values.
   *
   * \tparam Space: type of the discretization space on which the function
   * is defined, which must match the `LinearFunctionalSpaceConcept` concept
   * \tparam layout: data layout description, stating the number of
   * components of the function and the distance between the values of two
   * successive integration points, both of them being known at compile
   * time or only at run time
   * \tparam is_mutable: if true, the view allows modifying the values
   * stored in the external memory region, otherwise those values are constant
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout =
                FunctionDataLayoutDescription{},
            bool is_mutable = true>
  struct FunctionView
      : private PreconditionsChecker<FunctionView<Space, layout, is_mutable>>,
        public FunctionDataLayout<layout> {
    /*!
     * \brief type of the view to the external data used by the view, which
     * is a mutable span if `is_mutable` is true, a constant span otherwise
     */
    using ExternalData =
        std::conditional_t<is_mutable, std::span<real>, std::span<const real>>;
    /*!
     * \brief type of the value returned by the call operator, which is a
     * reference to a scalar value if the number of components is equal to
     * one at compile time, a mutable `std::span` otherwise
     */
    using ValuesView = std::conditional_t<
        layout.data_size == dynamic_extent,
        std::span<real>,
        std::conditional_t<layout.data_size == 1,
                           real&,
                           std::span<real, layout.data_size>>>;
    /*!
     * \brief type of the value returned by the call operator of a constant
     * view, i.e. a constant version of `ValuesView`
     */
    using ConstValuesView = std::conditional_t<
        layout.data_size == dynamic_extent,
        std::span<const real>,
        std::conditional_t<layout.data_size == 1,
                           const real&,
                           std::span<const real, layout.data_size>>>;
    /*!
     * \brief true if the view can be seen as a scalar valued function,
     * i.e. if the number of components is equal to one or is only known at
     * run time
     */
    static constexpr bool allowScalarAccessor =
        (layout.data_size == dynamic_extent) ? true : layout.data_size == 1;
    /*!
     * \brief check that the preconditions to build the view are met when
     * the number of components and the stride are known at compile time
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \return true if the external data is large enough to hold the values
     * of the function
     * \note a space of size zero is always accepted, which may happen due
     * to partitioning in parallel
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData) requires((layout.data_size != dynamic_extent) &&
                               (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met when
     * the number of components is only known at run time and the stride is
     * known at compile time
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     * \return true if the number of components is not greater than the
     * stride and the external data is large enough to hold the values of
     * the function
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met when
     * the number of components is known at compile time and the stride is
     * only known at run time
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dstride: stride of the data
     * \return true if the stride is not smaller than the number of
     * components and the external data is large enough to hold the values
     * of the function
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size != dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met when
     * both the number of components and the stride are only known at run
     * time, the stride being taken equal to the number of components
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     * \return true if the preconditions are met
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met when
     * both the number of components and the stride are only known at run time
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     * \param[in] dstride: stride of the data
     * \return true if the number of components and the stride are valid and
     * the external data is large enough to hold the values of the function
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met when
     * both the number of components and the stride are only known at run
     * time and are given by a data layout object
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout, which gives the number of components and
     * the stride of the function
     * \return true if the preconditions are met
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is known at compile
     * time and the stride is only known at run time
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dstride: stride of the data
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is known at compile
     * time and the stride is only known at run time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dstride: stride of the data
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     * \param[in] dstride: stride of the data
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     * \param[in] dstride: stride of the data
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time, the stride being taken equal to the
     * number of components
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time, the stride being taken equal to the
     * number of components
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is only known at run
     * time and the stride is known at compile time
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor when the number of components is only known at run
     * time and the stride is known at compile time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] dsize: number of components of the function
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are known at compile time
     * \param[in] s: quadrature space
     * \param[in] v: values
     */
    constexpr FunctionView(const Space&,
                           ExternalData)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are known at compile time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time and are given by a data layout object
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when both the number of components and the stride
     * are only known at run time and are given by a data layout object
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is known at compile
     * time and the stride is only known at run time and given by a data
     * layout object
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is known at compile
     * time and the stride is only known at run time and given by a data
     * layout object
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor when the number of components is only known at run
     * time and given by a data layout object, the stride being known at
     * compile time
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor when the number of components is only known at run
     * time and given by a data layout object, the stride being known at
     * compile time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: quadrature space
     * \param[in] v: values
     * \param[in] l: data layout
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief move constructor
     */
    constexpr FunctionView(FunctionView&&) = default;
    /*!
     * \brief copy constructor
     */
    constexpr FunctionView(const FunctionView&) = default;
    //! \return the underlying quadrature space
    [[nodiscard]] constexpr const Space& getSpace() const MGIS_NOEXCEPT;
    /*!
     * \brief return a pointer to the values associated with an integration
     * point, given its offset
     * \param[in] t: tag stating that no bound checking is performed
     * (`mgis::unsafe`), the offset must be strictly smaller than the size of
     * the space
     * \param[in] o: offset associated with the integration point
     * \return a pointer to the first value of the integration point
     * \note this function is available on mutable views defined on linear
     * element spaces without element workspaces
     */
    [[nodiscard]] constexpr real*
    data(mgis::attributes::UnsafeAttribute, const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \brief return a pointer to the values associated with an integration
     * point, given the index of its cell and its own index
     * \param[in] t: tag stating that no bound checking is performed
     * (`mgis::unsafe`), the indices must be valid
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return a pointer to the first value of the integration point
     * \note this function is available on mutable views defined on linear
     * quadrature spaces without cell workspaces
     */
    [[nodiscard]] constexpr real* data(mgis::attributes::UnsafeAttribute,
                                       const size_type,
                                       const size_type)  //
        requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \brief return a pointer to the constant values associated with an
     * integration point, given its offset
     * \param[in] t: tag stating that no bound checking is performed
     * (`mgis::unsafe`), the offset must be strictly smaller than the size of
     * the space
     * \param[in] o: offset associated with the integration point
     * \return a pointer to the first value of the integration point
     * \note this function is available on views defined on linear element
     * spaces without element workspaces
     */
    [[nodiscard]] constexpr const real* data(mgis::attributes::UnsafeAttribute,
                                             const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \brief return a pointer to the constant values associated with an
     * integration point, given the index of its cell and its own index
     * \param[in] t: tag stating that no bound checking is performed
     * (`mgis::unsafe`), the indices must be valid
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return a pointer to the first value of the integration point
     * \note this function is available on views defined on linear quadrature
     * spaces without cell workspaces
     */
    [[nodiscard]] constexpr const real* data(mgis::attributes::UnsafeAttribute,
                                             const size_type,
                                             const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \brief return the values associated with an integration point, given
     * its offset
     * \param[in] o: offset associated with the integration point
     * \return a view to the values of the integration point, of type
     * `ValuesView`
     * \note this function is available on mutable views defined on linear
     * element spaces without element workspaces
     */
    [[nodiscard]] constexpr ValuesView operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \brief return the values associated with an integration point, given
     * the index of its cell and its own index
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return a view to the values of the integration point, of type
     * `ValuesView`
     * \note this function is available on mutable views defined on linear
     * quadrature spaces without cell workspaces
     */
    [[nodiscard]] constexpr ValuesView
    operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \brief return the constant values associated with an integration
     * point, given its offset
     * \param[in] o: offset associated with the integration point
     * \return a view to the values of the integration point, of type
     * `ConstValuesView`
     * \note this function is available on views defined on linear element
     * spaces without element workspaces
     */
    [[nodiscard]] constexpr ConstValuesView operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \brief return the constant values associated with an integration
     * point, given the index of its cell and its own index
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return a view to the values of the integration point, of type
     * `ConstValuesView`
     * \note this function is available on views defined on linear quadrature
     * spaces without cell workspaces
     */
    [[nodiscard]] constexpr ConstValuesView operator()(const size_type,
                                                       const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \brief check if the given view can be used in place of the current one
     * \param[in] v: view to compare with the current view
     * \return true if the current view and the given view are defined on
     * equivalent spaces and have the same number of components
     */
    [[nodiscard]] constexpr bool checkCompatibility(const FunctionView&) const;
    /*!
     * \return a mutable view to all the values stored in the external
     * memory region
     * \note the returned span covers the whole external memory region, which
     * may contain data not related to the function, for instance when the
     * data stride is greater than the number of components
     */
    [[nodiscard]] constexpr std::span<real> data() requires(is_mutable);
    /*!
     * \return a view to all the constant values stored in the external
     * memory region
     * \note the returned span covers the whole external memory region, which
     * may contain data not related to the function, for instance when the
     * data stride is greater than the number of components
     */
    [[nodiscard]] constexpr std::span<const real> data() const;
    //! \brief destructor
    constexpr ~FunctionView() = default;

   protected:
    //! \brief underlying discretization space
    const Space space;
    //! \brief view to the values stored in the external memory region
    ExternalData values;
  };  // end of FunctionView

  /*!
   * \brief specialisation stating that a `FunctionView` is a lightweight
   * view, i.e. a function which is its own view: it can be copied at a
   * negligible cost
   *
   * \tparam Space: type of the discretization space
   * \tparam layout: data layout description
   * \tparam is_mutable: mutability flag of the view
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  struct LightweightViewTraits<FunctionView<Space, layout, is_mutable>>
      : std::true_type {};

  /*!
   * \brief an immutable `FunctionView`, i.e. a view which cannot modify the
   * values stored in the external memory region and which matches the
   * `EvaluatorConcept` concept
   *
   * \tparam Space: type of the discretization space
   * \tparam layout: data layout description
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout =
                FunctionDataLayoutDescription{}>
  using FunctionEvaluator = FunctionView<Space, layout, false>;

  /*!
   * \brief a noop function allowing `FunctionView` to match the
   * `EvaluatorConcept` concept
   * \tparam Space: type of the discretization space
   * \tparam layout: data layout description
   * \tparam is_mutable: mutability flag of the view
   * \param[in] eh: unused error handler
   * \param[in] e: view to be checked
   * \return always true: the values held by a view are always valid
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const FunctionView<Space, layout, is_mutable>&) MGIS_NOEXCEPT;
  /*!
   * \brief return the number of components of the given view
   * \tparam Space: type of the discretization space
   * \tparam layout: data layout description
   * \tparam is_mutable: mutability flag of the view
   * \param[in] e: view whose number of components is queried
   * \return the number of components of the values of the view
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr mgis::size_type getNumberOfComponents(
      const FunctionView<Space, layout, is_mutable>&) MGIS_NOEXCEPT;

  /*!
   * \brief an helper class storing the values of a function in a
   * `std::vector`
   *
   * \tparam Space: type of the discretization space on which the function
   * is defined
   * \tparam N: number of components, which may be `dynamic_extent`
   *
   * \note moving or copying the storage moves or copies the values but does
   * not update the views which point to them: this is the responsibility of
   * the classes deriving from this class
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  struct FunctionStorage {
    /*!
     * \brief constructor for a function whose number of components is known
     * at compile time: allocates `N * getSpaceSize(s)` values initialized
     * to zero
     * \param[in] s: space on which the function is defined
     */
    constexpr FunctionStorage(const Space& s) requires(N != dynamic_extent)
        : storage_values(N * getSpaceSize(s), real{}) {}
    /*!
     * \brief constructor for a function whose number of components is only
     * known at run time: allocates `getSpaceSize(s) * dsize` values
     * initialized to zero
     * \param[in] s: space on which the function is defined
     * \param[in] dsize: number of components of the function
     */
    constexpr FunctionStorage(const Space& s,
                              const size_type dsize) requires(N ==
                                                              dynamic_extent)
        : storage_values(getSpaceSize(s) * dsize, real{}) {}
    //! \brief move constructor
    constexpr FunctionStorage(FunctionStorage&&) = default;
    //! \brief copy constructor
    constexpr FunctionStorage(const FunctionStorage&) = default;
    //! \brief destructor
    constexpr ~FunctionStorage() = default;

   protected:
    //! \brief values
    std::vector<real> storage_values;
  };

  /*!
   * \brief data layout description whose number of components and stride
   * are both equal to `N`
   *
   * This variable is a workaround for a limitation of Visual Studio 2022 on
   * defining default values for non-type template arguments.
   *
   * \tparam N: number of components, which may be `dynamic_extent`
   */
  template <size_type N>
  inline constexpr auto simple_data_layout_description =
      FunctionDataLayoutDescription{.data_size = N, .data_stride = N};

  /*!
   * \brief a mutable `FunctionView` whose stride is equal to the number of
   * components, used to simplify the implementation of the `Function` class
   *
   * \tparam Space: type of the discretization space
   * \tparam N: number of components, which may be `dynamic_extent`
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  using FunctionViewAlias =
      FunctionView<Space, simple_data_layout_description<N>, true>;

  /*!
   * \brief default implementation of a function, which owns its values
   *
   * A `Function` stores its values in a `FunctionStorage` and derives from
   * the corresponding mutable `FunctionView`, so it can be called like a
   * function. It matches the `FunctionConcept` concept but not the
   * `EvaluatorConcept` concept, whose `check` function is deleted: a view
   * shall be created to evaluate a function.
   *
   * \tparam Space: type of the functional space
   * \tparam N: number of components, which may be `dynamic_extent` and
   * which must be greater than zero
   *
   * \note the data stride is equal to the data size
   */
  template <LinearFunctionalSpaceConcept Space, size_type N = dynamic_extent>
  requires(N > 0) struct Function
      : private PreconditionsChecker<Function<Space, N>>,
        public FunctionStorage<Space, N>,
        public FunctionViewAlias<Space, N> {
    /*!
     * \brief constructor from a space, for a function whose number of
     * components is known at compile time
     * \param[in] s: space on which the function is defined
     */
    constexpr Function(const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space and a data size, for a function whose
     * number of components is only known at run time
     * \param[in] s: space on which the function is defined
     * \param[in] dsize: number of components (data size)
     */
    constexpr Function(const Space&,
                       const size_type) requires(N == dynamic_extent);
    /*!
     * \brief check that the preconditions to build the function are met
     * when the number of components is known at compile time
     * \param[in] eh: unused error handler
     * \param[in] s: unused space
     * \return always true: no check is required in this case
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space, for a function whose number of
     * components is known at compile time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: space on which the function is defined
     */
    template <bool doPreconditionsCheck>
    constexpr Function(const PreconditionsCheck<doPreconditionsCheck>&,
                       const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief check that the preconditions to build the function are met
     * when the number of components is only known at run time
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] s: unused space
     * \param[in] dsize: number of components (data size)
     * \return true if the number of components is greater than zero
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        const size_type)  //
        requires(N == dynamic_extent);
    /*!
     * \brief constructor from a space and a data size, for a function whose
     * number of components is only known at run time
     * \tparam doPreconditionsCheck: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: space on which the function is defined
     * \param[in] dsize: data size
     */
    template <bool doPreconditionsCheck>
    constexpr Function(const PreconditionsCheck<doPreconditionsCheck>&,
                       const Space&,
                       const size_type) requires(N == dynamic_extent);
    /*!
     * \brief copy constructor when the number of components is only known
     * at run time: copies the values of the function and rebinds the view
     * to them
     * \param[in] f: function to copy
     */
    constexpr Function(const Function&) requires(N == dynamic_extent);
    /*!
     * \brief move constructor when the number of components is only known
     * at run time: moves the values of the function and rebinds the view to
     * them
     * \param[in] f: function to move
     */
    constexpr Function(Function&&) requires(N == dynamic_extent);
    /*!
     * \brief copy constructor when the number of components is known at
     * compile time: copies the values of the function and rebinds the view
     * to them
     * \param[in] f: function to copy
     */
    constexpr Function(const Function&) requires(N != dynamic_extent);
    /*!
     * \brief move constructor when the number of components is known at
     * compile time: moves the values of the function and rebinds the view
     * to them
     * \param[in] f: function to move
     */
    constexpr Function(Function&&) requires(N != dynamic_extent);
    /*!
     * \brief return a mutable view of the function
     * \return a mutable view sharing the values owned by the function
     */
    [[nodiscard]] constexpr FunctionView<Space,
                                         simple_data_layout_description<N>,
                                         true>
    view() &;
    /*!
     * \brief return an immutable view of the function
     * \return an immutable view of the function values
     */
    [[nodiscard]] constexpr FunctionView<Space,
                                         simple_data_layout_description<N>,
                                         false>
    view() const&;
    /*!
     * \brief deletion of the view of a rvalue function to prevent creating
     * a view pointing to values which are about to be destroyed
     */
    auto view() && = delete;
    //! \brief give access to the `data` functions of the base class
    using FunctionView<Space, simple_data_layout_description<N>, true>::data;
    /*!
     * \brief fill the structure using raw data
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] values: raw data, whose size must be equal to the total
     * number of values held by the function
     * \return true if the raw data has been copied in the function values
     */
    [[nodiscard]] constexpr bool fill(AbstractErrorHandler&,
                                      std::span<const real>) MGIS_NOEXCEPT;
    /*!
     * \brief fill the structure using raw data
     * \param[in] eh: error handler used to report the detected errors
     * \param[in] values: raw data, whose size must be equal to the total
     * number of values held by the function
     * \return true if the raw data has been copied in the function values
     */
    [[nodiscard]] constexpr bool fill(
        AbstractErrorHandler&, std::initializer_list<real>) MGIS_NOEXCEPT;
    //! \brief destructor
    constexpr ~Function();
  };

  /*!
   * \brief class template deduction guide: a `Function` constructed from a
   * space and a run-time number of components has a number of components
   * which is only known at run time
   * \tparam SpaceType: type of the discretization space
   */
  template <LinearFunctionalSpaceConcept SpaceType>
  Function(const SpaceType&, const size_type)
      -> Function<SpaceType, dynamic_extent>;

  /*!
   * \brief convert a function to a mutable view
   * \tparam Space: type of the functional space
   * \tparam N: number of components of the function
   * \param[in] f: function to convert
   * \return a mutable view of the function, which shares its values
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  [[nodiscard]] constexpr auto view(Function<Space, N>&);

  /*!
   * \brief convert a function to a mutable view
   * \tparam N: number of components of the view, which must be greater than
   * zero and known at compile time
   * \tparam Space: type of the functional space
   * \tparam N2: number of components of the function
   * \param[in] f: function to convert
   * \return a mutable view of the function, which shares its values
   * \note the `requires` clause makes this overload viable only when the
   * number of components of the function is known at compile time and
   * equal to `N`
   */
  template <size_type N, LinearFunctionalSpaceConcept Space, size_type N2>
  [[nodiscard]] constexpr auto view(Function<Space, N2>&)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2));

  /*!
   * \brief convert a function to an immutable view
   * \tparam Space: type of the functional space
   * \tparam N: number of components of the function
   * \param[in] f: function to convert
   * \return an immutable view of the function values
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  [[nodiscard]] constexpr auto view(const Function<Space, N>&);

  /*!
   * \brief deletion of the view of a rvalue function to prevent creating a
   * view pointing to values which are about to be destroyed
   * \tparam Space: type of the functional space
   * \tparam N: number of components of the function
   * \param[in] f: rvalue function
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  auto view(Function<Space, N>&&) = delete;

  /*!
   * \brief convert a function to an immutable view
   * \tparam N: number of components of the view, which must be greater than
   * zero and known at compile time
   * \tparam Space: type of the functional space
   * \tparam N2: number of components of the function
   * \param[in] f: function to convert
   * \return an immutable view of the function values
   * \note the `requires` clause makes this overload viable only when the
   * number of components of the function is known at compile time and
   * equal to `N`
   */
  template <size_type N, LinearFunctionalSpaceConcept Space, size_type N2>
  [[nodiscard]] constexpr auto view(const Function<Space, N2>&)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2));

  /*!
   * \brief deletion of the view of a rvalue function to prevent creating a
   * view pointing to values which are about to be destroyed
   * \tparam N: number of components of the view
   * \tparam Space: type of the functional space
   * \tparam N2: number of components of the function
   * \param[in] f: rvalue function
   */
  template <size_type N, LinearFunctionalSpaceConcept Space, size_type N2>
  auto view(Function<Space, N2>&&) = delete;

  /*!
   * \brief return the space on which the given view is defined
   * \tparam Space: type of the discretization space
   * \tparam layout: data layout description
   * \tparam is_mutable: mutability flag of the view
   * \param[in] f: view whose space is queried
   * \return the space on which the view is defined
   */
  template <LinearFunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr const auto& getSpace(
      const FunctionView<Space, layout, is_mutable>&);

  /*!
   * \brief deleted function so that `Function` does not match the
   * `EvaluatorConcept` concept: the values of a function must be evaluated
   * through a view
   * \tparam Space: type of the functional space
   * \tparam N: number of components of the function
   * \param[in] eh: error handler
   * \param[in] f: function to check
   */
  template <LinearFunctionalSpaceConcept Space, size_type N>
  bool check(AbstractErrorHandler&, const Function<Space, N>&) = delete;

}  // namespace mgis::function

#include "MGIS/Function/Function.ixx"

#endif /* LIB_MGIS_FUNCTION_FUNCTION_HXX */
