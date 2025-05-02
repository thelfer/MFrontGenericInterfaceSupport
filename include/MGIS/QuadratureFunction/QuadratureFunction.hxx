/*!
 * \file   MGIS/QuadratureFunction/QuadratureFunction.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_HXX

#include <span>
#include <limits>
#include <memory>
#include <vector>
#include <functional>
#include "MGIS/Config.hxx"

namespace mgis::quadrature_function {

  // forward declaration
  struct AbstractQuadratureSpace;

  /*!
   * \brief class describing the offset of the data manipulated by a quadrature
   * function.
   */
  template <size_type offset>
  struct QuadratureFunctionDataOffset {
    //! \return the offset of the first element
    constexpr size_type getDataOffset() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct QuadratureFunctionDataOffset<dynamic_extent> {
   protected:
    //! \return the offset of the first element
    size_type getDataOffset() const noexcept;

   protected:
    /*!
     * \brief begin of the data (offset of the function with respect to the
     * beginning of data values)
     */
    size_type data_begin = size_type{};
  };  // end of QuadratureFunctionDataOffset<dynamic_extent>

  /*!
   * \brief class describing the size (number of components) of the data
   * manipulated by a quadrature function.
   */
  template <size_type data_size>
  struct QuadratureFunctionDataSize {
    //! \return the number of components
    constexpr bool isScalar() const noexcept;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct QuadratureFunctionDataSize<dynamic_extent> {
    //! \return the number of components
    bool isScalar() const noexcept;
    //! \return the number of components
    size_type getNumberOfComponents() const noexcept;

   protected:
    //! \brief data size
    size_type data_size = size_type{};
  };  // end of QuadratureFunctionDataSize<dynamic_extent>

  /*!
   * \brief class describing the stride (number of components) of the data
   * manipulated by a quadrature function.
   */
  template <size_type data_stride>
  struct QuadratureFunctionDataStride {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    constexpr size_type getDataStride() const noexcept;
  };

  //! \brief partial specialisation in the dynamic_extent case
  template <>
  struct QuadratureFunctionDataStride<dynamic_extent> {
    /*!
     * \return the stride of data, i.e. the distance between the values of two
     * successive integration points.
     */
    size_type getDataStride() const noexcept;

   protected:
    //! \brief data size
    size_type data_stride = size_type{};
  };  // end of QuadratureFunctionDataStride<dynamic_extent>

  /*!
   * \brief a simple data structure describing how the data of a partial
   * quadrature function is mapped in memory
   */
  struct MGIS_EXPORT QuadratureFunctionDataLayout
      : public QuadratureFunctionDataOffset<dynamic_extent>,
        public QuadratureFunctionDataSize<dynamic_extent>,
        public QuadratureFunctionDataStride<dynamic_extent> {
    // \brief default constructor
    QuadratureFunctionDataLayout() = default;
    // \brief move constructor
    QuadratureFunctionDataLayout(QuadratureFunctionDataLayout&&) = default;
    // \brief copy constructor
    QuadratureFunctionDataLayout(const QuadratureFunctionDataLayout&) = default;
    // \brief move assignement
    QuadratureFunctionDataLayout& operator=(QuadratureFunctionDataLayout&&) =
        default;
    // \brief standard assignement
    QuadratureFunctionDataLayout& operator=(
        const QuadratureFunctionDataLayout&) = default;
    //
    using QuadratureFunctionDataOffset<dynamic_extent>::getDataOffset;
    //! \brief destructor
    ~QuadratureFunctionDataLayout() = default;

   protected:
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] i: integration point
     */
    size_type getDataOffset(const size_type) const noexcept;
  };  // end of struct QuadratureFunctionDataLayout

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   *
   * The `ImmutableQuadratureFunctionView` defines an immutable view
   * associated with a partial quadrature function on a memory region.
   *
   * This memory region may contain more data than the one associated with the
   * quadrature function as illustrated by the following figure:
   *
   * |---------------------------------------------------------------|
   * <-                         Raw data                            ->
   * |---------------------------------------|
   * <- Data of the first integration point-->
   * |      |---------------|                |
   *        <-function data->
   *        ^                                ^
   *        |                                |
   *    data_begin                           |
   *                                     data_size
   *
   * The size of the all data (including the one not related to the partial
   * quadrature function) associated with one integration point is called the
   * `data_stride` in the `ImmutableQuadratureFunctionView` class.
   *
   * Inside the data associated with on integration point, the function data
   * starts at the offset given by `data_begin`.
   *
   * The size of the data hold by the function per integration point, i.e. th
   * number of components of the function is given by `data_size`.
   */
  struct MGIS_EXPORT ImmutableQuadratureFunctionView
      : QuadratureFunctionDataLayout {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: offset of data
     * \param[in] ds: size of the data per integration points
     */
    ImmutableQuadratureFunctionView(
        std::shared_ptr<const AbstractQuadratureSpace>,
        std::span<const real>,
        const size_type = 0,
        const size_type = std::numeric_limits<size_type>::max());
    //! \return the underlying quadrature space
    const AbstractQuadratureSpace& getQuadratureSpace() const noexcept;
    //! \return the underlying quadrature space
    std::shared_ptr<const AbstractQuadratureSpace> getQuadratureSpacePointer()
        const noexcept;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */

    const real& getIntegrationPointValue(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type,
                                         const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    std::span<const real, N> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<const real> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<const real> getIntegrationPointValues(const size_type,
                                                    const size_type) const;
    /*!
     * \return if the current function has the same quadrature space and the
     * same number of components than the given view
     * \param[in] v: view
     */
    bool checkCompatibility(const ImmutableQuadratureFunctionView&) const;
    //! \return a view to the function values
    std::span<const real> getValues() const;
    //! \brief destructor
    ~ImmutableQuadratureFunctionView();

   protected:
    //! \brief default constructor
    ImmutableQuadratureFunctionView();
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] ds: data size
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view (stride)
     */
    ImmutableQuadratureFunctionView(
        std::shared_ptr<const AbstractQuadratureSpace>,
        const size_type,
        const size_type,
        const size_type);
    //! \brief underlying finite element space
    std::shared_ptr<const AbstractQuadratureSpace> qspace;
    //! \brief underlying values
    std::span<const real> immutable_values;
  };  // end of ImmutableQuadratureFunctionView

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   */
  struct MGIS_EXPORT QuadratureFunction : ImmutableQuadratureFunctionView {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] size: size of the data stored per integration points.
     */
    QuadratureFunction(std::shared_ptr<const AbstractQuadratureSpace>,
                       const size_type = 1);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view
     */
    QuadratureFunction(std::shared_ptr<const AbstractQuadratureSpace>,
                       std::span<real>,
                       const size_type = 0,
                       const size_type = std::numeric_limits<size_type>::max());
    /*!
     * \brief move constructor
     * \param[in] f: moved function
     * \param[in] local_copy: copy locally the function values if the moved
     * function does not holds them, i.e. is a view.
     * \note if the moved function holds the memory, the move constructor will
     * take ownership of the memory
     */
    QuadratureFunction(QuadratureFunction&&, const bool = false);
    //! \brief copy constructor
    QuadratureFunction(const QuadratureFunction&);
    /*!
     * \brief constructor from an immutable view
     *
     * \note: data are copied in a local array
     * \param[in] v: view
     */
    explicit QuadratureFunction(const ImmutableQuadratureFunctionView&);
    //! \brief assignement operator
    QuadratureFunction& operator=(const ImmutableQuadratureFunctionView&);
    //! \brief standard assignement operator
    QuadratureFunction& operator=(const QuadratureFunction&);
    //     //! \brief move assignement operator
    //     QuadratureFunction& operator=(QuadratureFunction&&);
    //
    using ImmutableQuadratureFunctionView::getIntegrationPointValue;
    using ImmutableQuadratureFunctionView::getIntegrationPointValues;
    using ImmutableQuadratureFunctionView::getValues;
    /*!
     * \brief return the value associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type);
    /*!
     * \brief return the value associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type, const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    std::span<real, N> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<real> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    std::span<real> getIntegrationPointValues(const size_type, const size_type);
    //! \brief destructor
    ~QuadratureFunction();

   protected:
    /*!
     * \brief turns this function into a view to the given function
     * \param[in] f: function
     */
    void makeView(QuadratureFunction&);
    /*!
     * \brief turns this function into a view to the given function
     * \param[in] f: function
     */
    void copy(const ImmutableQuadratureFunctionView&);
    /*!
     * \brief copy values for an immutable view
     * \param[in] v: view
     * \note: no compatibility checks are perfomed
     */
    void copyValues(const ImmutableQuadratureFunctionView&);
    //! \brief underlying values
    std::span<real> values;
    /*!
     * \brief storage for the values when the partial function holds the
     * values
     */
    std::vector<real> local_values_storage;
  };  // end of struct QuadratureFunction

}  // namespace mgis::quadrature_function

#include "MGIS/QuadratureFunction/QuadratureFunction.ixx"

#endif /* LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_HXX */
