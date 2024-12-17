/*!
 * \file   include/MGIS/MaterialProperty/OutputStatus.hxx
 * \brief
 * \author Thomas Helfer
 * \date   04/10/2022
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_MATERIALPROPERTY_OUTPUTSTATUS_HXX
#define LIB_MGIS_MATERIALPROPERTY_OUTPUTSTATUS_HXX

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
 * \brief this structure summarizes the exit status of a function conforming to
 * one of the `generic` material property interface.
 */
typedef struct {
  /*!
   * \brief exit status
   *
   * The exit status is zero if the result has been correctly evaluated.
   *
   * If the exit status is 1, a result has been computed, but it must be used
   * with caution. This is typically used to report that one argument was out of
   * its bounds.
   *
   * All negative values indicates that the result is not usable. For a material
   * property, the returned value is `nan`.
   *
   * For a material property, a negative value has the following meaning:
   *
   * - If the exit status is -1, an argument was out of its physical bounds, or
   *   out of its bounds and a strict out of bounds policy is declared.
   * - If the exit status is -2, a C++ exception was thrown. If the exception
   *   was a child of `std::exception`, the content of the string returned by
   *   the `what` method is copyied in the `message` field. Otherwise, the
   *   message field contains the "unknown exception" string.
   * - If the exit status is -3, an error occured in the `C` library, i.e. the
   *   `errno` value was set to a non zero value during the computation.
   *   The value of `errno` corresponding to the error is stored to in the
   *   `c_error_number` field of this structure. The string returned by
   *   `strerrno` is returned. Note that the `errno` value is always reset to
   *   the  value it had before the call.
   * - If the exit status is -4, the computed value is invalid (either \nan`,
   *   `inf`, or `-inf`).
   * - If the exit status is -5, the number of arguments is invalid.
   */
  int status;
  //! \brief error number reported by the C library.
  int c_error_number;
  /*!
   * \brief bounds status
   *
   * This status has the following meaning:
   *
   * - zero means that no argument was outside its bounds or its physical
   *   bounds.
   * - a negative values means that one argument went beyond its physical
   *   bounds. The absolute value gives the rank of this argument (here
   *   the rank starts at 1).
   * - a positive value means that one argument went beyond its bounds.
   *   The value gives the rank of this argument (here the rank starts at 1).
   */
  int bounds_status;
  //! \brief error message
  char msg[512];
} mgis_mp_OutputStatus;  // end of struct mgis_mp_OutputStatus

#ifdef __cplusplus
}  // end of extern "C"
#endif /* __cplusplus */

namespace mgis::material_property {

  //! \brief a simple alias
  using OutputStatus = ::mgis_mp_OutputStatus;

}  // end of namespace mgis::material_property

#endif /* LIB_MGIS_MATERIALPROPERTY_OUTPUTSTATUS_HXX */
