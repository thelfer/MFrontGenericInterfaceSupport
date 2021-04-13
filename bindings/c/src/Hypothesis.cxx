/*!
 * \file   bindings/c/src/Hypothesis.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/03/2020
 */

#include <cstring>
#include "MGIS/Behaviour/Hypothesis.h"

extern "C" {

mgis_status mgis_bv_get_space_dimension(mgis_size_type* const s,
                                        const char* const h) {
  if ((::strcmp(h, "AxisymmetricalGeneralisedPlaneStrain") == 0) ||
      (::strcmp(h, "AxisymmetricalGeneralisedPlaneStress") == 0)) {
    *s = 1u;
  } else if ((::strcmp(h, "Axisymmetrical") == 0) ||
             (::strcmp(h, "PlaneStress") == 0) ||
             (::strcmp(h, "PlaneStrain") == 0) ||
             (::strcmp(h, "GeneralisedPlaneStrain") == 0)) {
    *s = 2u;
  } else if (::strcmp(h, "Tridimensional") == 0) {
    *s = 3u;
  } else {
    return mgis_report_failure(
        "mgis_bv_get_space_dimension: "
        "unsupported hypothesis");
  }
  return mgis_report_success();
}  // end of mgis_bv_get_space_dimension

mgis_status mgis_bv_get_stensor_size(mgis_size_type* const s,
                                     const char* const h) {
  if ((::strcmp(h, "AxisymmetricalGeneralisedPlaneStrain") == 0) ||
      (::strcmp(h, "AxisymmetricalGeneralisedPlaneStress") == 0)) {
    *s = 3u;
  } else if ((::strcmp(h, "Axisymmetrical") == 0) ||
             (::strcmp(h, "PlaneStress") == 0) ||
             (::strcmp(h, "PlaneStrain") == 0) ||
             (::strcmp(h, "GeneralisedPlaneStrain") == 0)) {
    *s = 4u;
  } else if (::strcmp(h, "Tridimensional") == 0) {
    *s = 6u;
  } else {
    return mgis_report_failure(
        "mgis_bv_get_stensor_size: "
        "unsupported hypothesis");
  }
  return mgis_report_success();
}  // end of mgis_bv_get_stensor_size

mgis_status mgis_bv_get_tensor_size(mgis_size_type* const s,
                                    const char* const h) {
  if ((::strcmp(h, "AxisymmetricalGeneralisedPlaneStrain") == 0) ||
      (::strcmp(h, "AxisymmetricalGeneralisedPlaneStress") == 0)) {
    *s = 3u;
  } else if ((::strcmp(h, "Axisymmetrical") == 0) ||
             (::strcmp(h, "PlaneStress") == 0) ||
             (::strcmp(h, "PlaneStrain") == 0) ||
             (::strcmp(h, "GeneralisedPlaneStrain") == 0)) {
    *s = 5u;
  } else if (::strcmp(h, "Tridimensional") == 0) {
    *s = 9u;
  } else {
    return mgis_report_failure(
        "mgis_bv_get_tensor_size: "
        "unsupported hypothesis");
  }
  return mgis_report_success();
}  // end of mgis_bv_get_tensor_size

}  // end of extern "C"