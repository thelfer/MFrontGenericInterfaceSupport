/*!
 * \file   ThreadPool.h
 * \brief    
 * \author Thomas Helfer
 * \date   05/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_THREADPOOL_H
#define LIB_MGIS_THREADPOOL_H

#include "MGIS/Config.h"
#include "MGIS/Status.h"

#ifdef __cplusplus
#include "MGIS/ThreadPool.hxx"
#endif /*  __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

#ifdef __cplusplus
using mgis_ThreadPool = mgis::ThreadPool;
#else
/*!
 * \brief an opaque structure which can only be accessed through the MGIS' API.
 */
typedef struct mgis_ThreadPool mgis_ThreadPool;
#endif

/*!
 * \param[out] p: a pointer to the created thread pool
 * \param[in]  n: number of threads
 */
MGIS_C_EXPORT mgis_status mgis_create_thread_pool(mgis_ThreadPool**,
						  const mgis_size_type);
/*!
 * \param[in,out] p: a pointer to the thread pool to be destroyed
 */
MGIS_C_EXPORT mgis_status mgis_free_thread_pool(mgis_ThreadPool**);

#ifdef __cplusplus
} // end of extern "C"
#endif
  
#endif /* LIB_MGIS_THREADPOOL_H */
