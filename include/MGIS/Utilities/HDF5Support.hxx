/*!
 * \file   MGIS/Utilities/HDF5Support.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/01/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_UTILITIES_HDF5SUPPORT_HXX
#define LIB_MGIS_UTILITIES_HDF5SUPPORT_HXX

#include <span>
#include <string>
#include <vector>
#include <H5Cpp.h>
#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"

namespace mgis::utilities::hdf5 {

  template <typename T>
  H5::PredType getNativeType() noexcept;
  //
  template <>
  MGIS_EXPORT [[nodiscard]] H5::PredType getNativeType<float>() noexcept;
  template <>
  MGIS_EXPORT [[nodiscard]] H5::PredType getNativeType<double>() noexcept;
  template <>
  MGIS_EXPORT [[nodiscard]] H5::PredType getNativeType<long double>() noexcept;
  /*!
   * \brief check if an object with the given path exists in the given group
   * \param[in] g: group
   * \param[in] p: path
   */
  MGIS_EXPORT [[nodiscard]] bool exists(const H5::Group&,
                                        const std::string&) noexcept;
  /*!
   * \brief check if a group with the given path exists in the given given
   * \param[in] g: group
   * \param[in] p: path
   */
  MGIS_EXPORT [[nodiscard]] bool subGroupExists(const H5::Group&,
                                                const std::string&) noexcept;
  /*!
   * \brief create a new group
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group name
   */
  MGIS_EXPORT [[nodiscard]] std::optional<H5::Group> createGroup(
      Context&, const H5::Group&, const std::string&) noexcept;
  /*!
   * \brief open a new group
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group name
   */
  MGIS_EXPORT [[nodiscard]] std::optional<H5::Group> openGroup(
      Context&, const H5::Group&, const std::string&) noexcept;
  /*!
   * \brief open a data set
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: data set name
   */
  MGIS_EXPORT [[nodiscard]] std::optional<H5::DataSet> openDataSet(
      Context& ctx, const H5::Group&, const std::string&) noexcept;
  /*!
   * \brief remove a data set
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: data set name
   */
  MGIS_EXPORT [[nodiscard]] bool removeDataSet(Context&,
                                               const H5::Group&,
                                               const std::string&) noexcept;
  /*!
   * \param[in]  g: group
   * \param[in]  b: boolean allowing other objects than groups to be
   * inside the given group
   */
  MGIS_EXPORT [[nodiscard]] std::optional<std::vector<std::string>>
  getSubGroupNames(Context&, const H5::Group&, const bool) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] n: names
   * \param[in]  g: group
   * \param[in]  b: boolean allowing other objects than groups to be
   * inside the given group
   */
  MGIS_EXPORT [[nodiscard]] bool getSubGroupNames(Context&,
                                                  std::vector<std::string>&,
                                                  const H5::Group&,
                                                  const bool) noexcept;
  /*!
   * \return all the dataset names in a give group
   * \param[out, in] ctx: execution context
   * \param[in]  g: group
   */
  MGIS_EXPORT [[nodiscard]] std::optional<std::vector<std::string>>
  getDataSetNames(Context&, const H5::Group&) noexcept;
  /*!
   * \return all the dataset names in a give group
   * \param[out, in] ctx: execution context
   * \param[out] n: names
   * \param[in]  g: group
   */
  MGIS_EXPORT [[nodiscard]] bool getDataSetNames(Context&,
                                                 std::vector<std::string>&,
                                                 const H5::Group&) noexcept;
  /*!
   * \return true if the given group contains an object with the
   * given name
   * \param[out, in] ctx: execution context
   * \param[in] g: group
   * \param[in] n: name
   */
  MGIS_EXPORT [[nodiscard]] std::optional<bool> contains(
      Context&, const H5::Group&, const std::string&) noexcept;
  /*!
   * \brief delete an existing group or dataset if it exists
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group or dataset name
   */
  MGIS_EXPORT [[nodiscard]] bool unlinkIfExists(Context&,
                                                const H5::Group&,
                                                const std::string&) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param g : HDF5 group
   * \param n : name of the dataset
   * \param o : object to be written
   * \param[in]  b: allow overwrite
   */
  MGIS_EXPORT [[nodiscard]] bool write(Context&,
                                       H5::Group&,
                                       const std::string&,
                                       const real&,
                                       const bool) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be written
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT [[nodiscard]] bool read(Context&,
                                      real&,
                                      const H5::Group&,
                                      const std::string&) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param g : HDF5 group
   * \param n : name of the dataset
   * \param o : object to be written
   * \param[in]  b: allow overwrite
   */
  MGIS_EXPORT [[nodiscard]] bool write(Context&,
                                       H5::Group&,
                                       const std::string&,
                                       std::span<const real>,
                                       const bool) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be read
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT [[nodiscard]] bool read(Context&,
                                      std::vector<real>&,
                                      const H5::Group&,
                                      const std::string&) noexcept;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be read
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT [[nodiscard]] bool read(Context&,
                                      std::span<real>,
                                      const H5::Group&,
                                      const std::string&) noexcept;

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * a standard exception or an exception generated by the HD5 library
   *
   * \param[out] e: error back trace handler
   * \param[in] l: description of the call site
   */
  MGIS_EXPORT [[nodiscard]] InvalidResult registerH5ExceptionInErrorBacktrace(
      ErrorBacktrace&,
      const std::source_location& = std::source_location::current()) noexcept;
#else
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * a standard exception or an exception generated by the HD5 library
   *
   * \param[out] e: error back trace handler
   * \param[in] l: description of the call site
   */
  MGIS_EXPORT [[nodiscard]] InvalidResult registerH5ExceptionInErrorBacktrace(
      ErrorBacktrace&) noexcept;
#endif
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   *a standard exception or an exception generated by the HD5 library
   *
   *\param[out] e: error back trace handler
   */
  MGIS_EXPORT [[nodiscard]] InvalidResult
  registerH5ExceptionInErrorBacktraceWithoutSourceLocation(
      ErrorBacktrace&) noexcept;

}  // end of namespace mgis::utilities::hdf5

#endif /* LIB_MGIS_UTILITIES_HDF5SUPPORT_HXX */
