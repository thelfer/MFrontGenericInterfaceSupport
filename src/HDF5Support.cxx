/*!
 * \file   src/HDF5Support.cxx
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

#include <utility>
#include "MGIS/Utilities/HDF5Support.hxx"

namespace mgis::utilities::hdf5 {

  InvalidResult registerH5ExceptionInErrorBacktraceWithoutSourceLocation(
      ErrorBacktrace& e) noexcept {
    try {
      throw;
    } catch (H5::GroupIException& exception) {
      std::ignore =
          e.registerErrorMessageWithoutSourceLocation(exception.getDetailMsg());
    } catch (H5::Exception& exception) {
      std::ignore =
          e.registerErrorMessageWithoutSourceLocation(exception.getDetailMsg());
    } catch (std::exception& exception) {
      std::ignore = e.registerErrorMessageWithoutSourceLocation(
          std::string{exception.what()});
    } catch (...) {
      std::ignore = e.registerErrorMessageWithoutSourceLocation(
          "unknown exception thrown");
    }
    return {};
  }

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  InvalidResult registerH5ExceptionInErrorBacktrace(
      ErrorBacktrace& e, const std::source_location& l) noexcept {
    try {
      throw;
    } catch (H5::GroupIException& exception) {
      std::ignore = e.registerErrorMessage(exception.getDetailMsg());
    } catch (H5::Exception& exception) {
      std::ignore = e.registerErrorMessage(exception.getDetailMsg());
    } catch (std::exception& exception) {
      std::ignore = e.registerErrorMessage(std::string{exception.what()}, l);
    } catch (...) {
      std::ignore = e.registerErrorMessage("unknown exception thrown");
    }
    return {};
  }
#else
  InvalidResult registerH5ExceptionInErrorBacktrace(
      ErrorBacktrace& e) noexcept {
    return registerH5ExceptionInErrorBacktraceWithoutSourceLocation(e);
  }
#endif

  static std::vector<std::string> tokenize(std::string_view s,
                                           const char c,
                                           const bool keep_empty_strings) {
    std::vector<std::string> res;
    auto b = std::string::size_type{};
    auto e = s.find_first_of(c, b);
    while (std::string::npos != e || std::string::npos != b) {
      // Found a token, add it to the vector.
      res.push_back(std::string{s.substr(b, e - b)});
      if (keep_empty_strings) {
        b = e == std::string::npos ? e : e + 1;
      } else {
        b = s.find_first_not_of(c, e);
      }
      e = s.find_first_of(c, b);
    }
    return res;
  }  // end of tokenize

  bool exists(const H5::Group& g, const std::string& p) noexcept {
    // break paths
    const auto leafs = tokenize(p, '/', false);
    if (leafs.empty()) {
      return false;
    }
    auto cpath = std::string{};
    for (const auto& l : leafs) {
      if (!cpath.empty()) {
        cpath += '/';
      }
      cpath += l;
      const auto b = H5Lexists(g.getId(), cpath.c_str(), H5P_DEFAULT) > 0;
      if (!b) {
        return false;
      }
    }
    return true;
  }

  bool subGroupExists(const H5::Group& g, const std::string& p) noexcept {
    if (!exists(g, p)) {
      return false;
    }
    H5O_info_t infobuf;
#ifdef H5Oget_info_by_name_vers
#if H5Oget_info_by_name_vers >= 3
    auto status = H5Oget_info_by_name(g.getId(), p.c_str(), &infobuf,
                                      H5O_INFO_ALL, H5P_DEFAULT);
#else
    auto status =
        H5Oget_info_by_name(g.getId(), p.c_str(), &infobuf, H5P_DEFAULT);
#endif
#else
    auto status =
        H5Oget_info_by_name(g.getId(), p.c_str(), &infobuf, H5P_DEFAULT);
#endif
    return (status >= 0) && (infobuf.type == H5O_TYPE_GROUP);
  }

  std::optional<H5::Group> createGroup(Context& ctx,
                                       const H5::Group& g,
                                       const std::string& n) noexcept {
    if (subGroupExists(g, n)) {
      return openGroup(ctx, g, n);
    }
    try {
      H5::Group gr = g.createGroup(n);
      return gr;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of createGroup

  std::optional<H5::Group> openGroup(Context& ctx,
                                     const H5::Group& g,
                                     const std::string& n) noexcept {
    try {
      H5::Group gr = g.openGroup(n);
      return gr;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of openGroup

  bool removeDataSet(Context& ctx,
                     const H5::Group& g,
                     const std::string& n) noexcept {
    try {
      g.unlink(n);
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of removeDataSet

  std::optional<H5::DataSet> openDataSet(Context& ctx,
                                         const H5::Group& g,
                                         const std::string& n) noexcept {
    try {
      H5::DataSet d = g.openDataSet(n);
      return d;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of openDataSet

  std::optional<std::vector<std::string>> getSubGroupNames(
      Context& ctx, const H5::Group& g, const bool b) noexcept {
    std::vector<std::string> names;
    if (!getSubGroupNames(ctx, names, g, b)) {
      return {};
    }
    return names;
  }  // end of getSubGroupNames

  bool getSubGroupNames(Context& ctx,
                        std::vector<std::string>& n,
                        const H5::Group& g,
                        const bool b) noexcept {
    n.clear();
    try {
      const hsize_t s = g.getNumObjs();
      for (hsize_t i = 0; i != s; ++i) {
        if (g.getObjTypeByIdx(i) == H5G_GROUP) {
          n.push_back(g.getObjnameByIdx(i));
        } else if (!b) {
          return ctx.registerErrorMessage("getSubGroupNames: object '" +
                                          g.getObjnameByIdx(i) +
                                          "' is not a group");
        }
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of getSubGroupNames

  std::optional<std::vector<std::string>> getDataSetNames(
      Context& ctx, const H5::Group& g) noexcept {
    std::vector<std::string> names;
    if (!getDataSetNames(ctx, names, g)) {
      return {};
    }
    return names;
  }

  bool getDataSetNames(Context& ctx,
                       std::vector<std::string>& n,
                       const H5::Group& g) noexcept {
    n.clear();
    try {
      const hsize_t s = g.getNumObjs();
      for (hsize_t i = 0; i != s; ++i) {
        if (g.getObjTypeByIdx(i) == H5G_DATASET) {
          n.push_back(g.getObjnameByIdx(i));
        }
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of getDataSetNames

  std::optional<bool> contains(Context& ctx,
                               const H5::Group& g,
                               const std::string& n) noexcept {
    try {
      bool found = false;
      const hsize_t s = g.getNumObjs();
      for (hsize_t i = 0; (i != s) && (!found); ++i) {
        if (g.getObjnameByIdx(i) == n) {
          found = true;
        }
      }
      return found;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of contains

  bool unlinkIfExists(Context& ctx,
                      const H5::Group& g,
                      const std::string& n) noexcept {
    try {
      if (exists(g, n)) {
        g.unlink(n);
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }

  template <>
  [[nodiscard]] H5::PredType getNativeType<float>() noexcept {
    return H5::PredType::NATIVE_FLOAT;
  }  // end of getNativeType

  template <>
  [[nodiscard]] H5::PredType getNativeType<double>() noexcept {
    return H5::PredType::NATIVE_DOUBLE;
  }  // end of getNativeType

  template <>
  [[nodiscard]] H5::PredType getNativeType<long double>() noexcept {
    return H5::PredType::NATIVE_LDOUBLE;
  }  // end of getNativeType

  MGIS_EXPORT [[nodiscard]] bool write(Context& ctx,
                                       H5::Group& g,
                                       const std::string& d,
                                       const real& o,
                                       const bool b) noexcept {
    if (b) {
      if (!unlinkIfExists(ctx, g, d)) {
        return false;
      }
    }
    try {
      hsize_t dimsf[1] = {1};
      H5::DataSpace dataspace(1, dimsf);
      auto dataset = g.createDataSet(d, getNativeType<real>(), dataspace);
      dataset.write(&o, getNativeType<real>());
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of write

  MGIS_EXPORT [[nodiscard]] bool read(Context& ctx,
                                      real& o,
                                      const H5::Group& g,
                                      const std::string& d) noexcept {
    auto odataset = openDataSet(ctx, g, d);
    if (isInvalid(odataset)) {
      return false;
    }
    try {
      const auto s = odataset->getSpace();
      if (s.getSimpleExtentNdims() != 1) {
        return ctx.registerErrorMessage("madnex::read: invalid type size");
      }
      hsize_t dims[1];
      s.getSimpleExtentDims(dims);
      if (dims[0] != 1) {
        return ctx.registerErrorMessage("madnex::read: invalid type");
      }
      odataset->read(&o, getNativeType<real>());
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of read

  bool write(Context& ctx,
             H5::Group& g,
             const std::string& d,
             std::span<const real> o,
             const bool b) noexcept {
    if (b) {
      if (!unlinkIfExists(ctx, g, d)) {
        return false;
      }
    }
    try {
      if (o.empty()) {
        auto c = real{};
        hsize_t dimsf[1];
        dimsf[0] = 1u;
        H5::DataSpace dataspace(1, dimsf);
        H5::StrType datatype(getNativeType<real>());
        auto dataset = g.createDataSet(d, datatype, dataspace);
        dataset.write(&c, getNativeType<real>());
      } else {
        hsize_t dimsf[1];
        dimsf[0] = o.size();
        H5::DataSpace dataspace(1, dimsf);
        auto dataset = g.createDataSet(d, getNativeType<real>(), dataspace);
        dataset.write(&o[0], getNativeType<real>());
      }
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }

  bool read(Context& ctx,
            std::vector<real>& o,
            const H5::Group& g,
            const std::string& d) noexcept {
    auto odataset = openDataSet(ctx, g, d);
    if (isInvalid(odataset)) {
      return false;
    }
    try {
      H5::DataSpace filespace = odataset->getSpace();
      hsize_t dims[1];
      filespace.getSimpleExtentDims(dims);
      o.resize(dims[0]);
      odataset->read(o.data(), getNativeType<real>());
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }

  bool read(Context& ctx,
            std::span<real> values,
            const H5::Group& g,
            const std::string& n) noexcept {
    using namespace mgis::utilities::hdf5;
    const auto odataset = openDataSet(ctx, g, n);
    if (isInvalid(odataset)) {
      return false;
    }
    try {
      H5::DataSpace filespace = odataset->getSpace();
      hsize_t dims[1];
      filespace.getSimpleExtentDims(dims);
      if (values.size() != dims[0]) {
        return ctx.registerErrorMessage("invalid size while reading dataset '" +
                                        n + "'");
      }
      odataset->read(values.data(), getNativeType<real>());
      return true;
    } catch (...) {
      std::ignore = registerH5ExceptionInErrorBacktrace(ctx);
    }
    return false;
  }  // end of restore

}  // end of namespace mgis::utilities::hdf5