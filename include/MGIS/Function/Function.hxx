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
#include <memory>
#include <vector>
#include "MGIS/Config.hxx"
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
    constexpr bool isScalar() const noexcept;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct FunctionDataSize<dynamic_extent> {
    //! \return the number of components
    bool isScalar() const noexcept;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;

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
    constexpr size_type getDataStride() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct FunctionDataStride<dynamic_extent> {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    size_type getDataStride() const noexcept;

   protected:
    //! \brief data size
    size_type data_stride = size_type{};
  };  // end of FunctionDataStride<dynamic_extent>

  struct DataLayoutDescription {
    size_type data_size = dynamic_extent;
    size_type data_stride = dynamic_extent;
  };

  //! \brief a simple helper function
  constexpr bool has_dynamic_properties(const DataLayoutDescription&);

  /*!
   * \brief a simple data structure describing how the data of a
   * function is mapped in memory
   */
  template <DataLayoutDescription layout>
  requires((layout.data_size > 0) &&
           (layout.data_stride > 0)) struct MGIS_EXPORT DataLayout
      : public FunctionDataSize<layout.data_size>,
        public FunctionDataStride<layout.data_stride> {
    // \brief default constructor
    DataLayout() = default;
    // \brief move constructor
    DataLayout(DataLayout&&) = default;
    // \brief copy constructor
    DataLayout(const DataLayout&) = default;
    // \brief move assignement
    DataLayout& operator=(DataLayout&&) = default;
    // \brief standard assignement
    DataLayout& operator=(const DataLayout&) = default;
    //! \brief destructor
    ~DataLayout() = default;

   protected:
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] i: integration point
     */
    size_type getDataOffset(const size_type) const noexcept;
  };  // end of struct DataLayout

  /*!
   * \brief function defined on a space.
   *
   * The `FunctionView ` defines a function using an external memory region.
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
   * the `FunctionView ` class.
   *
   * The size of the data hold by the function per element of the space, i.e. th
   * number of components of the function is given by `data_size`.
   */
  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout = {},
            bool is_mutable = true>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      struct FunctionView : DataLayout<layout> {
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
     * \param[in] ctx: execution context.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        std::shared_ptr<const Space>,
        ExternalData) requires((layout.data_size != dynamic_extent) &&
                               (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] ctx: execution context.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        std::shared_ptr<const Space>,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] ctx: execution context.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dstride: data stride
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        std::shared_ptr<const Space>,
        ExternalData,
        const size_type) requires((layout.data_size != dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] ctx: execution context.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        std::shared_ptr<const Space>,
        ExternalData,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] ctx: execution context.
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dsize: size of the data per elements
     * \param[in] dstride: data stride
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        std::shared_ptr<const Space>,
        ExternalData,
        const size_type,
        const size_type) requires((layout.data_size == dynamic_extent) &&
                                  (layout.data_stride == dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] dstride: size of the data per elements
     */
    FunctionView(std::shared_ptr<const Space>,
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
    FunctionView(std::shared_ptr<const Space>,
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
    FunctionView(std::shared_ptr<const Space>,
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
    FunctionView(std::shared_ptr<const Space>,
                 ExternalData,
                 const size_type)  //
        requires((layout.data_size == dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    FunctionView(std::shared_ptr<const Space>,
                 ExternalData)  //
        requires((layout.data_size != dynamic_extent) &&
                 (layout.data_stride != dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    FunctionView(std::shared_ptr<const Space>,
                 ExternalData,
                 const DataLayout<layout>&) requires((layout.data_size ==
                                                      dynamic_extent) &&
                                                     (layout.data_stride ==
                                                      dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    FunctionView(std::shared_ptr<const Space>,
                 ExternalData,
                 const DataLayout<layout>&) requires((layout.data_size !=
                                                      dynamic_extent) &&
                                                     (layout.data_stride ==
                                                      dynamic_extent));
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] l: data layout
     */
    FunctionView(std::shared_ptr<const Space>,
                 ExternalData,
                 const DataLayout<layout>&) requires((layout.data_size ==
                                                      dynamic_extent) &&
                                                     (layout.data_stride !=
                                                      dynamic_extent));
    //! \return the underlying quadrature space
    const Space& getSpace() const noexcept;
    //! \return the underlying quadrature space
    std::shared_ptr<const Space> getSpacePointer() const noexcept;
    //! \brief a noop function to match the EvaluatorConcept concept
    bool check(Context&) const noexcept;
    //! \brief a noop function to match the EvaluatorConcept concept
    void allocateWorkspace() noexcept;
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    real* data(mgis::attributes::UnsafeAttribute, const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    real* data(mgis::attributes::UnsafeAttribute,
               const size_type,
               const size_type)  //
        requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    const real* data(mgis::attributes::UnsafeAttribute, const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    const real* data(mgis::attributes::UnsafeAttribute,
                     const size_type,
                     const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    ValuesView operator()(const size_type) requires(
        is_mutable&& LinearElementSpaceConcept<Space> &&
        (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    ValuesView operator()(const size_type, const size_type) requires(
        is_mutable&& LinearQuadratureSpaceConcept<Space> &&
        (!hasCellWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    ConstValuesView operator()(const size_type) const
        requires(LinearElementSpaceConcept<Space> &&
                 (!hasElementWorkspace<Space>));
    /*!
     * \return the data associated with an integration point
     * \param[in] e: element index
     * \param[in] i: quadrature point index
     */
    ConstValuesView operator()(const size_type, const size_type) const
        requires(LinearQuadratureSpaceConcept<Space> &&
                 (!hasCellWorkspace<Space>));
    /*!
     * \return if the current function has the same quadrature space and the
     * same number of components than the given view
     * \param[in] v: view
     */
    bool checkCompatibility(const FunctionView&) const;
    //! \return a view to the function values
    std::span<const real> data() const;
    //! \brief destructor
    ~FunctionView() noexcept;

   protected:
    //! \brief underlying discretization space
    std::shared_ptr<const Space> space;
    //! \brief underlying values
    ExternalData values;
  };  // end of FunctionView

  //! \brief a simple alias
  template <FunctionalSpaceConcept Space, DataLayoutDescription layout = {}>
  using FunctionEvaluator = FunctionView<Space, layout, false>;

  /*!
   * \brief an helper class storing the values of a function
   */
  template <FunctionalSpaceConcept Space, size_type N>
  struct FunctionStorage {
    FunctionStorage(const Space& s) requires(N != dynamic_extent)
        : storage_values(N * s.size()) {}
    FunctionStorage(const Space& s,
                    const size_type dsize) requires(N == dynamic_extent)
        : storage_values(s.size() * dsize) {}
    FunctionStorage(FunctionStorage&&) = default;
    FunctionStorage(const FunctionStorage&) = default;

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
      : public FunctionStorage<Space, N>,
        public FunctionView<Space, {.data_size = N, .data_stride = N}, true> {
    /*!
     * \brief constructor from a space and a data size
     * \param[in] ctx: execution context
     * \param[in] s: space
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        const std::shared_ptr<const Space>&)  //
        requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space
     * \param[in] s: space
     */
    Function(std::shared_ptr<const Space>) requires(N != dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] ctx: execution context
     * \param[in] s: space
     * \param[in] dsize: data size
     */
    [[nodiscard]] static bool checkPreconditions(
        Context&,
        const std::shared_ptr<const Space>&,
        const size_type)  //
        requires(N == dynamic_extent);
    /*!
     * \brief constructor from a space and a data size
     * \param[in] s: space
     * \param[in] dsize: data size
     */
    Function(std::shared_ptr<const Space>,
             const size_type) requires(N == dynamic_extent);
    //! \brief copy constructor
    Function(const Function&) requires(N == dynamic_extent);
    //! \brief assignement constructor
    Function(Function&&) requires(N == dynamic_extent);
    //! \brief copy constructor
    Function(const Function&) requires(N != dynamic_extent);
    //! \brief assignement constructor
    Function(Function&&) requires(N != dynamic_extent);
    //! \brief return a view of the function
    FunctionView<Space, {.data_size = N, .data_stride = N}, true> view();
    //! \brief return a view of the function
    FunctionView<Space, {.data_size = N, .data_stride = N}, false> view() const;

   protected:
    /*
     * This function is made protected to avoid Function from being treated
     * as an evaluator
     */
    using FunctionView<Space, {.data_size = N, .data_stride = N}, true>::check;
    /*
     * This function is made protected to avoid Function from being treated as
     * an evaluator
     */
    using FunctionView<Space, {.data_size = N, .data_stride = N}, true>::
        allocateWorkspace;
  };

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <FunctionalSpaceConcept Space, size_type N>
  auto view(const Function<Space, N>&);

  /*!
   * \brief convert a function to a immutable view
   * \param[in] f: function
   */
  template <size_type N, FunctionalSpaceConcept Space, size_type N2>
  auto view(const Function<Space, N2>&)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2));

}  // namespace mgis::function

#include "MGIS/Function/Function.ixx"

#endif /* LIB_MGIS_FUNCTION_FUNCTION_HXX */
