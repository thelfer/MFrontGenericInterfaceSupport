/*!
 * \file   Variable.h
 * \brief    
 * \author Thomas Helfer
 * \date   25/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_VARIABLE_H
#define LIB_MGIS_BEHAVIOUR_VARIABLE_H

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

/*!
 * \brief type of a variable
 */
typedef enum {
  MGIS_BV_SCALAR  = 0,
  MGIS_BV_VECTOR  = 1,
  MGIS_BV_STENSOR = 2,
  MGIS_BV_TENSOR  = 3
} mgis_bv_VariableType;

/*!
 * \return the size of an array able to store the variable
 * \param [out] s: size
 * \param [in] h: modelling hypothesis
 * \param [in] t: variable type
 */
MGIS_C_EXPORT mgis_status mgis_bv_get_variable_size(mgis_size_type* const,
						    const char* const,
						    const mgis_bv_VariableType);

#ifdef __cplusplus
}  // end of extern "C"
#endif

#endif /* LIB_MGIS_BEHAVIOUR_VARIABLE_H */
