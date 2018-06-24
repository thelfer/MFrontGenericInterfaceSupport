/*!
 * \file   Status.h
 * \brief
 * \author Thomas Helfer
 * \date   24/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_STATUS_H
#define LIB_MGIS_STATUS_H

#include "MGIS/Config-c.h"

#ifdef __cplusplus
extern "C" {
#endif /*  __cplusplus */

typedef enum { MGIS_SUCCESS, MGIS_FAILURE } mgis_exit_flag;

/*!
 * structure in charge of reporting if a function call has been successfull or
 * not. In case of failure, the `msg` fields can be used to retrieve an error
 * message.
 */
typedef struct {
  //! exit status
  enum { SUCCESS, FAILURE } exit_status;
  /*!
   * \brief a pointer to a per thread buffer in which the error message is
   * stored.
   * This buffer is allocated statically and shall not be freed by the caller.
   */
  const char *msg;
} mgis_status; // end of struct mgis_status

/*!
 * \brief build a success status
 * \return the build status
 */
MGIS_C_EXPORT mgis_status mgis_report_success();

/*!
 * \brief build a failed status
 * \param[in] e: error message
 * \note the error message is copied in a statically allocated buffer.
 * \return the build status
 */
MGIS_C_EXPORT mgis_status mgis_report_failure(const char *const);

/*!
 * \brief build a status from the current c++ exception
 * \note if the exception is derived from `std::exception`, the error message
 * returned by the `what` method is copied in a statically allocated buffer.
 * \note this function is meant to be used internally and is not exported.
 * \return the build status
 */
mgis_status mgis_handle_cxx_exception();

#ifdef __cplusplus
}
#endif /*  __cplusplus */

#endif /* LIB_MGIS_STATUS_H */
