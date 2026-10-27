/*!
 * \file   bindings/c/src/Variable.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/03/2020
 */

#include "MGIS/Behaviour/Hypothesis.h"
#include "MGIS/Behaviour/Variable.h"

extern "C" {

mgis_status mgis_bv_get_variable_size(mgis_size_type* const s,
                                      const char* const h,
                                      const mgis_bv_VariableType t) {
  if (t == MGIS_BV_SCALAR) {
    *s = 1u;
    return mgis_report_success();
  } else if (t == MGIS_BV_VECTOR) {
    return mgis_bv_get_space_dimension(s, h);
  } else if (t == MGIS_BV_VECTOR_1D) {
    *s = 1u;
    return mgis_report_success();
  } else if (t == MGIS_BV_VECTOR_2D) {
    *s = 2u;
    return mgis_report_success();
  } else if (t == MGIS_BV_VECTOR_3D) {
    *s = 3u;
    return mgis_report_success();
  } else if (t == MGIS_BV_STENSOR) {
    return mgis_bv_get_stensor_size(s, h);
  } else if (t == MGIS_BV_STENSOR_1D) {
    *s = 3u;
    return mgis_report_success();
  } else if (t == MGIS_BV_STENSOR_2D) {
    *s = 4u;
    return mgis_report_success();
  } else if (t == MGIS_BV_STENSOR_3D) {
    *s = 6u;
    return mgis_report_success();
  } else if (t == MGIS_BV_TENSOR) {
    return mgis_bv_get_tensor_size(s, h);
  } else if (t == MGIS_BV_TENSOR_1D) {
    *s = 3u;
    return mgis_report_success();
  } else if (t == MGIS_BV_TENSOR_2D) {
    *s = 5u;
    return mgis_report_success();
  } else if (t == MGIS_BV_TENSOR_3D) {
    *s = 9u;
    return mgis_report_success();
  }
  return mgis_report_failure(
      "mgis_bv_get_variable_size: "
      "unsupported variable type");
}  // end of mgis_bv_get_variable_size

}  // end of extern "C"
