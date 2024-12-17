/*!
 * \file   include/MGIS/MaterialProperty/OutOfBoundsPolicy.hxx
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

#ifndef LIB_MGIS_MATERIALPROPERTY_OUTOFBOUNDSPOLICY_HXX
#define LIB_MGIS_MATERIALPROPERTY_OUTOFBOUNDSPOLICY_HXX

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*!
 * \brief available out of bounds policies
 */
typedef enum {
  MGIS_MATERIALPROPERTY_NONE_POLICY,    /*!<
                                         * With this policy, nothing is done if
                                         * the arguments are    out of their
                                         * bounds    (checks are not even
                                         * performed).
                                         */
  MGIS_MATERIALPROPERTY_WARNING_POLICY, /*!<
                                         * With this policy, checks on the
                                         * arguments are performed. If one
                                         * argument if out of its bounds,
                                         * this will be reported in the
                                         * output status and an
                                         * appropriate error message will be
                                         * reported. The computations are
                                         * however performed.
                                         */
  MGIS_MATERIALPROPERTY_STRICT_POLICY   /*!<
                                         * With this policy, checks on the
                                         * arguments are   performed. If one
                                         * argument   if out of its bounds,
                                         * this   will be reported in the
                                         * output   status and an   appropriate
                                         * error   message will be reported.
                                         */
} mgis_mp_OutOfBoundsPolicy;            // end of mfront_gmp_OutOfBoundsPolicy

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus

namespace mgis::material_property {

  using OutOfBoundsPolicy = ::mgis_mp_OutOfBoundsPolicy;

}  // end of namespace mgis::material_property

#endif /* __cplusplus */

#endif /* LIB_MGIS_MATERIALPROPERTY_OUTOFBOUNDSPOLICY_HXX */
