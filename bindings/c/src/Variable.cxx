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
  } else if (t == MGIS_BV_STENSOR) {
    return mgis_bv_get_stensor_size(s, h);
  } else if (t == MGIS_BV_TENSOR) {
    return mgis_bv_get_tensor_size(s, h);
  }
  return mgis_report_failure(
      "mgis_bv_get_variable_size: "
      "unsupported variable type");
}  // end of mgis_bv_get_variable_size

}  // end of extern "C"
