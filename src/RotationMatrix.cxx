/*!
 * \file   src/RotationMatrix.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   17/02/2021
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/RotationMatrix.hxx"

namespace mgis::behaviour {

  static std::vector<mgis::real> copyValuesIfRequired(
      const mgis::span<const mgis::real>& values, const mgis::StorageMode& s) {
    if (s != mgis::StorageMode::LOCAL_STORAGE) {
      return {};
    }
    return {values.begin(), values.end()};
  }

  static mgis::span<const mgis::real> initializeLocalSpan(
      const mgis::span<const mgis::real>& evalues,
      const std::vector<mgis::real>& lvalues,
      const mgis::StorageMode& s) {
    if (s == mgis::StorageMode::LOCAL_STORAGE) {
      return lvalues;
    }
    return evalues;
  }  // end of initializeLocalSpan

  MaterialAxisStorage::MaterialAxisStorage(
      const mgis::span<const mgis::real>& values, const mgis::StorageMode& s)
      : a_values(copyValuesIfRequired(values, s)),
        a(initializeLocalSpan(values, this->a_values, s)) {
  }  // end of MaterialAxisStorage::MaterialAxisStorage

  MaterialAxisStorage::~MaterialAxisStorage() = default;

  static mgis::span<const mgis::real> checkMaterialAxis2D(
      const mgis::span<const mgis::real>& v) {
    if (v.empty()) {
      mgis::raise(
          "RotationMatrix2D::RotationMatrix2D: "
          "empty values for material axis in 2D");
    }
    if (!(v.size() % 2)) {
      mgis::raise(
          "RotationMatrix2D::RotationMatrix2D: "
          "invalid number of values for material axis in 2D");
    }
    return v;
  } // end of checkMaterialAxis3D

  RotationMatrix2D::RotationMatrix2D(const mgis::span<const mgis::real>& v,
                                     const mgis::StorageMode& s)
      : MaterialAxisStorage(checkMaterialAxis2D(v), s) {}

  RotationMatrix2D::~RotationMatrix2D() = default;

  static mgis::span<const mgis::real> checkMaterialAxis3D(
      const mgis::span<const mgis::real>& v) {
    if (v.empty()) {
      mgis::raise(
          "RotationMatrix3D::RotationMatrix3D: "
          "empty values for material axis in 3D");
    }
    const auto s = v.size();
    if (!(s % 3)) {
      mgis::raise(
          "RotationMatrix3D::RotationMatrix3D: "
          "invalid number of values for material axis in 3D");
    }
    return v;
  } // end of checkMaterialAxis3D

  RotationMatrix3D::RotationMatrix3D(const mgis::span<const mgis::real>& v1,
                                     const mgis::StorageMode& s1,
                                     const mgis::span<const mgis::real>& v2,
                                     const mgis::StorageMode& s2)
      : a1(checkMaterialAxis3D(v1), s1), a2(checkMaterialAxis3D(v2), s2) {
    if (v1.empty()) {
      mgis::raise(
          "RotationMatrix3D::RotationMatrix3D: "
          "empty values for the first axis");
    }
  }

  RotationMatrix3D::RotationMatrix3D(const mgis::span<const mgis::real>& v1,
                                     const mgis::span<const mgis::real>& v2,
                                     const mgis::StorageMode& s)
      : a1(checkMaterialAxis3D(v1), s), a2(checkMaterialAxis3D(v2), s) {}

  RotationMatrix3D::~RotationMatrix3D() = default;

} // end of namespace mgis::behaviour