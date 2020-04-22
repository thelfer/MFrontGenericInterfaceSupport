/*!
 * \file   Hypothesis.h
 * \brief    
 * \author Thomas Helfer
 * \date   02/03/2020
 */

#ifndef LIB_MGIS_BEHAVIOUR_HYPOTHESIS_H
#define LIB_MGIS_BEHAVIOUR_HYPOTHESIS_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

/*!
 * \brief get the size of a tensor for the given modelling hypothesis
 * \param [out] s: size
 * \param [in] h: modelling hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_space_dimension(mgis_size_type* const,
                                                      const char* const);
/*!
 * \brief get the size of a tensor for the given modelling hypothesis
 * \param [out] s: size
 * \param [in] h: modelling hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_stensor_size(mgis_size_type* const,
                                                  const char*const);
/*!
 * \brief get the size of a tensor for the given modelling hypothesis
 * \param [out] s: size
 * \param [in] h: modelling hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_tensor_size(mgis_size_type* const,
                                                  const char* const);

#ifdef __cplusplus
}  // end of extern "C"
#endif

#endif /* LIB_MGIS_BEHAVIOUR_HYPOTHESIS_H */
