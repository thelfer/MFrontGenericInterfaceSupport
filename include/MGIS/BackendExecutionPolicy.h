/*!
 * \file   include/MGIS/BackendExecutionPolicy.h
 * \brief  This file declares the type of the execution policies of the
 * backends of the `generic-parallel` interfaces of `MFront`. Those policies
 * are defined by `TFEL` (`MFront/GenericParallel/` headers, base class
 * `mfront::generic_parallel::AbstractBackendExecutionPolicy`): `MGIS` only
 * gives them to the generated functions, so that a declaration is enough and
 * `MGIS` does not depend on `TFEL` at compile time. In `C`, the policies are
 * opaque and can only be given as null pointers (default policy of the
 * backend).
 * \date   22/10/2026
 * \copyright Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BACKENDEXECUTIONPOLICY_H
#define LIB_MGIS_BACKENDEXECUTIONPOLICY_H

// NOLINTBEGIN(modernize-use-using)
#ifdef __cplusplus

namespace mfront::generic_parallel {

  //! \brief base class of the execution policies of the backends
  struct AbstractBackendExecutionPolicy;

}  // end of namespace mfront::generic_parallel

//! \brief a simple alias
typedef mfront::generic_parallel::AbstractBackendExecutionPolicy
    mgis_gp_AbstractBackendExecutionPolicy;

#else /* __cplusplus */

//! \brief opaque structure
typedef struct mgis_gp_AbstractBackendExecutionPolicy
    mgis_gp_AbstractBackendExecutionPolicy;

#endif /* __cplusplus */
// NOLINTEND(modernize-use-using)

#endif /* LIB_MGIS_BACKENDEXECUTIONPOLICY_H */
