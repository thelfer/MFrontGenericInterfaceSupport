/*!
 * \file   bindings/c/include/MGIS/Model.h
 * \brief    
 * \author Thomas Helfer
 * \date   14/10/2021
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_MODEL_MODEL_H
#define LIB_MGIS_MODEL_MODEL_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"
#include "MGIS/Behaviour/Behaviour.h"

#ifdef __cplusplus
#include "MGIS/Model/Model.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
//! \brief a simple alias
using mgis_model_Model = mgis::model::Model;
#else
//! \brief a simple alias
typedef mgis_bv_Behaviour mgis_model_Model;
#endif

/*!
 * \brief load a behaviour
 *
 * \param[out] ptr: behaviour
 * \param[in] l: library name
 * \param[in] m: model name
 * \param[in] h: hypothesis
 */
MGIS_C_EXPORT mgis_status mgis_model_load(mgis_model_Model**,
                                          const char* const,
                                          const char* const,
                                          const char* const);
/*!
 * \brief free the memory associated with the given model.
 * \param[in,out] m: model
 */
MGIS_C_EXPORT mgis_status mgis_model_free_model(mgis_model_Model**);

#ifdef __cplusplus
} // end of extern "C"
#endif /*  __cplusplus */

#endif /* LIB_MGIS_MODEL_MODEL_H */
