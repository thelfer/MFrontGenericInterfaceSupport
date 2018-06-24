/*!
 * \file   Config-c.h
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

#ifndef LIB_CONFIG_C_H
#define LIB_CONFIG_C_H

#include "MGIS/Config.h"

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericInterface_c_EXPORTS
#define MGIS_C_EXPORT MGIS_VISIBILITY_EXPORT
#else
#ifndef MGIS_STATIC_BUILD
#define MGIS_C_EXPORT MGIS_VISIBILITY_IMPORT
#else
#define MGIS_C_EXPORT
#endif
#endif
#else
#define MGIS_C_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

#endif /* LIB_CONFIG_C_H */
