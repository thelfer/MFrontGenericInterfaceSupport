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
#include "MGIS/Function/Space.hxx"

namespace mgis::function {

  /*!
   * \brief class describing the size (number of components) of the data
   * manipulated by a quadrature function.
   */
  template <size_type data_size>
  struct FunctionDataSize {
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
  struct FunctionDataStride {
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
    size_type size = dynamic_extent;
    size_type stride = dynamic_extent;
  };

  //! \brief a simple helper function
  constexpr bool has_dynamic_properties(const DataLayoutDescription&);

  /*!
   * \brief a simple data structure describing how the data of a partial
   * quadrature function is mapped in memory
   */
  template <DataLayoutDescription layout>
  struct MGIS_EXPORT DataLayout
      : public FunctionDataSize<layout.size>,
        public FunctionDataStride<layout.stride> {
    //
    static_assert((layout.size > 0) && (layout.stride >= layout.size),
                  "invalid layout");
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
   * \brief quadrature function defined on a partial quadrature space.
   *
   * The `ImmutableFunctionView` defines an immutable view
   * associated with a partial quadrature function on a memory region.
   *
   * This memory region may contain more data than the one associated with the
   * quadrature function as illustrated by the following figure:
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
   * The size of the all data (including the one not related to the partial
   * quadrature function) associated with one integration point is called the
   * `data_stride` in the `ImmutableFunctionView` class.
   *
   * Inside the data associated with on integration point, the function data
   * starts at the offset given by `data_begin`.
   *
   * The size of the data hold by the function per integration point, i.e. th
   * number of components of the function is given by `data_size`.
   */
  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  struct MGIS_EXPORT ImmutableFunctionView : DataLayout<layout> {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] ds: size of the data per elements
     */
    ImmutableFunctionView(std::shared_ptr<const Space>,
                          std::span<const real>,
                          const size_type)  //
        requires(layout.size == dynamic_extent);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     */
    ImmutableFunctionView(std::shared_ptr<const Space>,
                          std::span<const real>);
    //! \return the underlying quadrature space
    const Space& getSpace() const noexcept;
    //! \return the underlying quadrature space
    std::shared_ptr<const Space> getSpacePointer() const noexcept;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getValue(const size_type) const
        requires(LinearSpaceConcept<Space>);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    std::span<const real, N> getValues(const size_type) const
        requires(LinearSpaceConcept<Space>);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    std::span<const real> getValues(const size_type) const
        requires(LinearSpaceConcept<Space>);
    /*!
     * \return if the current function has the same quadrature space and the
     * same number of components than the given view
     * \param[in] v: view
     */
    bool checkCompatibility(const ImmutableFunctionView&) const;
    //! \return a view to the function values
    std::span<const real> getValues() const;
    //! \brief destructor
    ~ImmutableFunctionView();

   protected:
    //! \brief default constructor
    ImmutableFunctionView();
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] ds: data size
     * \param[in] ds: size of the view (stride)
     */
    ImmutableFunctionView(
        std::shared_ptr<const Space>,
        const size_type,
        const size_type) requires((layout.size != dynamic_extent) &&
                                  (layout.stride != dynamic_extent));
    //! \brief underlying discretization space
    std::shared_ptr<const Space> space;
    //! \brief underlying values
    std::span<const real> immutable_values;
  };  // end of ImmutableFunctionView

  //   /*!
  //    * \brief quadrature function defined on a partial quadrature space.
  //    */
  //   struct MGIS_EXPORT Function : ImmutableFunctionView {
  //     /*!
  //      * \brief constructor
  //      * \param[in] s: quadrature space.
  //      * \param[in] size: size of the data stored per integration points.
  //      */
  //     Function(std::shared_ptr<const Space>, const size_type = 1);
  //     /*!
  //      * \brief constructor
  //      * \param[in] s: quadrature space.
  //      * \param[in] v: values
  //      * \param[in] ds: size of the view
  //      */
  //     Function(std::shared_ptr<const Space>,
  //              std::span<real>,
  //              const size_type = std::numeric_limits<size_type>::max());
  //     /*!
  //      * \brief move constructor
  //      * \param[in] f: moved function
  //      * \param[in] local_copy: copy locally the function values if the moved
  //      * function does not holds them, i.e. is a view.
  //      * \note if the moved function holds the memory, the move constructor
  //      will
  //      * take ownership of the memory
  //      */
  //     Function(Function&&, const bool = false);
  //     //! \brief copy constructor
  //     Function(const Function&);
  //     /*!
  //      * \brief constructor from an immutable view
  //      *
  //      * \note: data are copied in a local array
  //      * \param[in] v: view
  //      */
  //     explicit Function(const ImmutableFunctionView&);
  //     //! \brief assignement operator
  //     Function& operator=(const ImmutableFunctionView&);
  //     //! \brief standard assignement operator
  //     Function& operator=(const Function&);
  //     //     //! \brief move assignement operator
  //     //     Function& operator=(Function&&);
  //     //
  //     using ImmutableFunctionView::getValue;
  //     using ImmutableFunctionView::getValues;
  //     /*!
  //      * \brief return the value associated with an integration point
  //      * \param[in] o: offset associated with the integration point
  //      * \note this method is only meaningful when the quadrature function is
  //      * scalar
  //      */
  //     real& getValue(const size_type);
  //     /*!
  //      * \brief return the value associated with an integration point
  //      * \param[in] e: global element number
  //      * \param[in] i: integration point number in the element
  //      * \note this method is only meaningful when the quadrature function is
  //      * scalar
  //      */
  //     real& getValue(const size_type, const size_type);
  //     /*!
  //      * \brief return the data associated with an integration point
  //      * \param[in] o: offset associated with the integration point
  //      */
  //     template <size_type N>
  //     std::span<real, N> getValues(const size_type);
  //     /*!
  //      * \brief return the data associated with an integration point
  //      * \param[in] o: offset associated with the integration point
  //      */
  //     std::span<real> getValues(const size_type);
  //     //! \brief destructor
  //     ~Function();
  //
  //    protected:
  //     /*!
  //      * \brief turns this function into a view to the given function
  //      * \param[in] f: function
  //      */
  //     void makeView(Function&);
  //     /*!
  //      * \brief turns this function into a view to the given function
  //      * \param[in] f: function
  //      */
  //     void copy(const ImmutableFunctionView&);
  //     /*!
  //      * \brief copy values for an immutable view
  //      * \param[in] v: view
  //      * \note: no compatibility checks are perfomed
  //      */
  //     void copyValues(const ImmutableFunctionView&);
  //     //! \brief underlying values
  //     std::span<real> values;
  //     /*!
  //      * \brief storage for the values when the partial function holds the
  //      * values
  //      */
  //     std::vector<real> local_values_storage;
  //   };  // end of struct Function

}  // namespace mgis::function

#include "MGIS/Function/Function.ixx"

#endif /* LIB_MGIS_FUNCTION_FUNCTION_HXX */
