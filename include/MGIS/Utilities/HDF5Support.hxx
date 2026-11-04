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
#include "MGIS/Utilities/HDF5Forward.hxx"

namespace mgis::utilities::hdf5 {

  template <typename T>
  H5::PredType getNativeType() MGIS_NOEXCEPT;
  //
  template <>
  MGIS_EXPORT MGIS_NODISCARD H5::PredType getNativeType<float>() MGIS_NOEXCEPT;
  template <>
  MGIS_EXPORT MGIS_NODISCARD H5::PredType getNativeType<double>() MGIS_NOEXCEPT;
  template <>
  MGIS_EXPORT MGIS_NODISCARD H5::PredType getNativeType<long double>()
      MGIS_NOEXCEPT;
  /*!
   * \brief check if an object with the given path exists in the given group
   * \param[in] g: group
   * \param[in] p: path
   */
  MGIS_EXPORT MGIS_NODISCARD bool exists(const H5::Group&,
                                         const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief check if a group with the given path exists in the given given
   * \param[in] g: group
   * \param[in] p: path
   */
  MGIS_EXPORT MGIS_NODISCARD bool subGroupExists(
      const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief create a new group
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group name
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<H5::Group> createGroup(
      Context&, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief open a new group
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group name
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<H5::Group> openGroup(
      Context&, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief open a data set
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: data set name
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<H5::DataSet> openDataSet(
      Context& ctx, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief remove a data set
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: data set name
   */
  MGIS_EXPORT MGIS_NODISCARD bool removeDataSet(
      Context&, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \param[in]  g: group
   * \param[in]  b: boolean allowing other objects than groups to be
   * inside the given group
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<std::vector<std::string>>
  getSubGroupNames(Context&, const H5::Group&, const bool) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] n: names
   * \param[in]  g: group
   * \param[in]  b: boolean allowing other objects than groups to be
   * inside the given group
   */
  MGIS_EXPORT MGIS_NODISCARD bool getSubGroupNames(Context&,
                                                   std::vector<std::string>&,
                                                   const H5::Group&,
                                                   const bool) MGIS_NOEXCEPT;
  /*!
   * \return all the dataset names in a give group
   * \param[out, in] ctx: execution context
   * \param[in]  g: group
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<std::vector<std::string>>
  getDataSetNames(Context&, const H5::Group&) MGIS_NOEXCEPT;
  /*!
   * \return true if the given group contains an object with the
   * given name
   * \param[out, in] ctx: execution context
   * \param[in] g: group
   * \param[in] n: name
   */
  MGIS_EXPORT MGIS_NODISCARD std::optional<bool> contains(
      Context&, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \brief delete an existing group or dataset if it exists
   * \param[out, in] ctx: execution context
   * \param[in] g: parent group
   * \param[in] n: group or dataset name
   */
  MGIS_EXPORT MGIS_NODISCARD bool unlinkIfExists(
      Context&, const H5::Group&, const std::string&) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param g : HDF5 group
   * \param n : name of the dataset
   * \param o : object to be written
   * \param[in]  b: allow overwrite
   */
  MGIS_EXPORT MGIS_NODISCARD bool write(Context&,
                                        H5::Group&,
                                        const std::string&,
                                        const real&,
                                        const bool) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be written
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT MGIS_NODISCARD bool read(Context&,
                                       real&,
                                       const H5::Group&,
                                       const std::string&) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param g : HDF5 group
   * \param n : name of the dataset
   * \param o : object to be written
   * \param[in]  b: allow overwrite
   */
  MGIS_EXPORT MGIS_NODISCARD bool write(Context&,
                                        H5::Group&,
                                        const std::string&,
                                        std::span<const real>,
                                        const bool) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be read
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT MGIS_NODISCARD bool read(Context&,
                                       std::vector<real>&,
                                       const H5::Group&,
                                       const std::string&) MGIS_NOEXCEPT;
  /*!
   * \param[out, in] ctx: execution context
   * \param[out] o : object to be read
   * \param g : HDF5 group
   * \param n : name of the dataset
   */
  MGIS_EXPORT MGIS_NODISCARD bool read(Context&,
                                       std::span<real>,
                                       const H5::Group&,
                                       const std::string&) MGIS_NOEXCEPT;

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * a standard exception or an exception generated by the HD5 library
   *
   * \param[out] e: error back trace handler
   * \param[in] l: description of the call site
   */
  MGIS_EXPORT MGIS_NODISCARD InvalidResult registerH5ExceptionInErrorBacktrace(
      ErrorBacktrace&,
      const std::source_location& = std::source_location::current())
      MGIS_NOEXCEPT;
#else
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * a standard exception or an exception generated by the HD5 library
   *
   * \param[out] e: error back trace handler
   * \param[in] l: description of the call site
   */
  MGIS_EXPORT MGIS_NODISCARD InvalidResult
  registerH5ExceptionInErrorBacktrace(ErrorBacktrace&) MGIS_NOEXCEPT;
#endif
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   *a standard exception or an exception generated by the HD5 library
   *
   *\param[out] e: error back trace handler
   */
  MGIS_EXPORT MGIS_NODISCARD InvalidResult
  registerH5ExceptionInErrorBacktraceWithoutSourceLocation(ErrorBacktrace&)
      MGIS_NOEXCEPT;

}  // end of namespace mgis::utilities::hdf5

#endif /* LIB_MGIS_UTILITIES_HDF5SUPPORT_HXX */
