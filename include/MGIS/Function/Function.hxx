/*!
 * \file   MGIS/Function/Function.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
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
   * manipulated by a quadrature function.
   */
  template <size_type data_size>
  requires(data_size > 0) struct FunctionDataSize {
    //! \return the number of components
    [[nodiscard]] constexpr bool isScalar() const noexcept;
    //! \return the number of components
    [[nodiscard]] constexpr size_type getNumberOfComponents() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct FunctionDataSize<dynamic_extent> {
    //! \return the number of components
    [[nodiscard]] constexpr bool isScalar() const noexcept;
    //! \return the number of components
    [[nodiscard]] constexpr size_type getNumberOfComponents() const noexcept;

   protected:
    //! \brief data size
    size_type data_size = size_type{};
  };  // end of FunctionDataSize<dynamic_extent>

  /*!
   * \brief class describing the stride (number of components) of the data
   * manipulated by a quadrature function.
   */
  template <size_type data_stride>
  requires(data_stride > 0) struct FunctionDataStride {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    [[nodiscard]] constexpr size_type getDataStride() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct FunctionDataStride<dynamic_extent> {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    [[nodiscard]] constexpr size_type getDataStride() const noexcept;

   protected:
    //! \brief data size
    size_type data_stride = size_type{};
  };  // end of FunctionDataStride<dynamic_extent>

  struct FunctionDataLayoutDescription {
    size_type data_size = dynamic_extent;
    size_type data_stride = dynamic_extent;
  };

  //! \brief a simple helper function
  constexpr bool has_dynamic_properties(const FunctionDataLayoutDescription&);

  /*!
   * \brief a simple data structure describing how the data of a
   * function is mapped in memory
   */
  template <FunctionDataLayoutDescription layout>
  requires((layout.data_size > 0) &&
           (layout.data_stride > 0)) struct MGIS_EXPORT FunctionDataLayout
      : public FunctionDataSize<layout.data_size>,
        public FunctionDataStride<layout.data_stride> {
    // \brief default constructor
    FunctionDataLayout() = default;
    // \brief move constructor
    FunctionDataLayout(FunctionDataLayout&&) = default;
    // \brief copy constructor
    FunctionDataLayout(const FunctionDataLayout&) = default;
    // \brief move assignement
    FunctionDataLayout& operator=(FunctionDataLayout&&) = default;
    // \brief standard assignement
    FunctionDataLayout& operator=(const FunctionDataLayout&) = default;
    //! \brief destructor
    ~FunctionDataLayout() = default;

   protected:
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] i: integration point
     */
    constexpr size_type getDataOffset(const size_type) const noexcept;
  };  // end of struct FunctionDataLayout

  /*!
   * \brief function defined on a space.
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
   * |---------------|                       |
   * <-function data->
   *                 ^                       ^
   *                 |                       |
   *             data_size                   |
   *                                     data_stride
   *
   * The size of the all data (including the one not related to the function)
   * associated with one element of the space is called the `data_stride` in
   * the `FunctionView` class.
   *
   * The size of the data hold by the function per element of the space, i.e. th
   * number of components of the function is given by `data_size`.
   */
  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout = {},
            bool is_mutable = true>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      struct FunctionView
      : private PreconditionsChecker<FunctionView<Space, layout, is_mutable>>,
        public FunctionDataLayout<layout> {
    //! \brief a simple alias to the type holding the data used by the view
    using ExternalData =
        std::conditional_t<is_mutable, std::span<real>, std::span<const real>>;
    //! \brief type of the value returned by the call operator
    using ValuesView = std::conditional_t<
        layout.data_size == dynamic_extent,
        std::span<real>,
        std::conditional_t<layout.data_size == 1,
                           real&,
                           std::span<real, layout.data_size>>>;
    //! \brief type of the value returned by the call operator (const case)
    using ConstValuesView = std::conditional_t<
        layout.data_size == dynamic_extent,
        std::span<const real>,
        std::conditional_t<layout.data_size == 1,
                           const real&,
                           std::span<const real, layout.data_size>>>;
    //
    static constexpr bool allowScalarAccessor =
        (layout.data_size == dynamic_extent) ? true : layout.data_size == 1;
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData) requires((layout.data_size != dynamic_extent) &&
                               (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dstride: data stride
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size != dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     * \param[in] dstride: data stride
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const size_type,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     * \param[in] dstride: data stride
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        ExternalData,
        const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dstride: size of the data per elements
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dstride: size of the data per elements
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     * \param[in] dstride: data stride
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     * \param[in] dstride: data stride
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
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData,
                           const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    constexpr FunctionView(const Space&,
                           ExternalData)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    template <bool doPreconditionsCheck>
    constexpr FunctionView(const PreconditionsCheck<doPreconditionsCheck>&,
                           const Space&,
                           ExternalData)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
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
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
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
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    constexpr FunctionView(const Space&,
                           ExternalData,
                           const FunctionDataLayout<layout>&)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
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
    //
    constexpr FunctionView(FunctionView&&) = default;
    constexpr FunctionView(const FunctionView&) = default;
    //! \return the underlying quadrature space
    [[nodiscard]] constexpr const Space& getSpace() const noexcept;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr real*
    data(mgis::attributes::UnsafeAttribute, const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr real* data(mgis::attributes::UnsafeAttribute,
                                       const size_type,
                                       const size_type)  //
        requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr const real* data(mgis::attributes::UnsafeAttribute,
                                             const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr const real* data(mgis::attributes::UnsafeAttribute,
                                             const size_type,
                                             const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr ValuesView operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ValuesView
    operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    [[nodiscard]] constexpr ConstValuesView operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    [[nodiscard]] constexpr ConstValuesView operator()(const size_type,
                                                       const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return if the current function has the same quadrature space and the
     * same number of components than the given view
     * \param[in] v: view
     */
    [[nodiscard]] constexpr bool checkCompatibility(const FunctionView&) const;
    //! \return a view to the function values
    [[nodiscard]] constexpr std::span<const real> data() const;
    //! \brief destructor
    constexpr ~FunctionView() = default;

   protected:
    //! \brief underlying discretization space
    const Space space;
    //! \brief underlying values
    ExternalData values;
  };  // end of FunctionView

  //! \brief a simple alias
  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout = {}>
  using FunctionEvaluator = FunctionView<Space, layout, false>;

  //! \brief a noop function to match the EvaluatorConcept concept
  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const FunctionView<Space, layout, is_mutable>&) noexcept;
  /*!
   * \brief do nothing function to match the Evaluator concept
   * \param[in] v: view
   */
  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  constexpr void allocateWorkspace(
      FunctionView<Space, layout, is_mutable>&) noexcept;
  /*!
   * \return the number of components of a view
   * \param[in] v: view
   */
  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr mgis::size_type getNumberOfComponents(
      const FunctionView<Space, layout, is_mutable>&) noexcept;

  /*!
   * \brief an helper class storing the values of a function
   */
  template <FunctionalSpaceConcept Space, size_type N>
  struct FunctionStorage {
    constexpr FunctionStorage(const Space& s) requires(N != dynamic_extent)
        : storage_values(N * getSpaceSize(s), real{}) {}
    constexpr FunctionStorage(const Space& s,
                              const size_type dsize) requires(N ==
                                                              dynamic_extent)
        : storage_values(getSpaceSize(s) * dsize, real{}) {}
    constexpr FunctionStorage(FunctionStorage&&) = default;
    constexpr FunctionStorage(const FunctionStorage&) = default;
    constexpr ~FunctionStorage() = default;

   protected:
    //! \brief values
    std::vector<real> storage_values;
  };

  /*!
   * \brief default implementation of a function
   *
   * \tparam Space: type of the functional space
   * \tparam N: number of components
   *
   * \note the data stride is equal to the data size
   */
  template <FunctionalSpaceConcept Space, size_type N = dynamic_extent>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      struct Function
      : private PreconditionsChecker<Function<Space, N>>,
        public FunctionStorage<Space, N>,
        public FunctionView<Space, {.data_size = N, .data_stride = N}, true> {
    /*!
     * \brief constructor from a space
     * \param[in] s: space
     */
    constexpr Function(const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] s: space
     * \param[in] dsize: number of components (data size)
     */
    constexpr Function(const Space&,
                       const size_type) requires(N == dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] eh: error handler
     * \param[in] s: space
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space
     * \param[in] s: space
     */
    template <bool doPreconditionsCheck>
    constexpr Function(const PreconditionsCheck<doPreconditionsCheck>&,
                       const Space&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] eh: error handler
     * \param[in] s: space
     * \param[in] dsize: data size
     */
    [[nodiscard]] static constexpr bool checkPreconditions(
        AbstractErrorHandler&,
        const Space&,
        const size_type)  //
        requires(N == dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] s: space
     * \param[in] dsize: data size
     */
    template <bool doPreconditionsCheck>
    constexpr Function(const PreconditionsCheck<doPreconditionsCheck>&,
                       const Space&,
                       const size_type) requires(N == dynamic_extent);
    //! \brief copy constructor
    constexpr Function(const Function&) requires(N == dynamic_extent);
    //! \brief assignement constructor
    constexpr Function(Function&&) requires(N == dynamic_extent);
    //! \brief copy constructor
    constexpr Function(const Function&) requires(N != dynamic_extent);
    //! \brief assignement constructor
    constexpr Function(Function&&) requires(N != dynamic_extent);
    //! \brief return a view of the function
    [[nodiscard]] constexpr FunctionView<Space,
                                         {.data_size = N, .data_stride = N},
                                         true>
    view();
    //! \brief return a view of the function
    [[nodiscard]] constexpr FunctionView<Space,
                                         {.data_size = N, .data_stride = N},
                                         false>
    view() const;
    //
    using FunctionView<Space, {.data_size = N, .data_stride = N}, true>::data;
    //! \return a view to the function values
    [[nodiscard]] constexpr std::span<real> data();
    /*!
     * \brief fill the structure using raw data
     *
     * \param[in] eh: error handler
     * \param[in] values: raw data
     */
    [[nodiscard]] constexpr bool fill(AbstractErrorHandler&,
                                      std::span<const real>) noexcept;
    /*!
     * \brief fill the structure using raw data
     *
     * \param[in] eh: error handler
     * \param[in] values: raw data
     */
    [[nodiscard]] constexpr bool fill(AbstractErrorHandler&,
                                      std::initializer_list<real>) noexcept;
    //! \brief destructor
    constexpr ~Function();
  };

  // class template deduction guide
  template <FunctionalSpaceConcept SpaceType>
  Function(const SpaceType&, const size_type)
      -> Function<SpaceType, dynamic_extent>;

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <FunctionalSpaceConcept Space, size_type N>
  [[nodiscard]] constexpr auto view(const Function<Space, N>&);

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space, size_type N2>
  [[nodiscard]] constexpr auto view(const Function<Space, N2>&)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2));

  template <FunctionalSpaceConcept Space,
            FunctionDataLayoutDescription layout,
            bool is_mutable>
  [[nodiscard]] constexpr const auto& getSpace(
      const FunctionView<Space, layout, is_mutable>&);

  //! \brief deleted function so that Function does not match the
  //! EvaluatorConcept
  template <FunctionalSpaceConcept Space, size_type N>
  bool check(AbstractErrorHandler&, const Function<Space, N>&) = delete;
  //! \brief deleted function so that Function does not match the
  //! EvaluatorConcept
  template <FunctionalSpaceConcept Space, size_type N>
  void allocateWorkspace(const Function<Space, N>&) = delete;

}  // namespace mgis::function

#include "MGIS/Function/Function.ixx"

#endif /* LIB_MGIS_FUNCTION_FUNCTION_HXX */
