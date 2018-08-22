/*!
 * \file   include/MGIS/Config.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_CONFIG_HXX
#define LIB_MGIS_CONFIG_HXX

#include"MGIS/Config-c.h"

namespace mgis {

  //! a simple alias
  using size_type = mgis_size_type;

  //! alias to the numeric type used
  using real = mgis_real;

}  // end of namespace mgis

#endif /* LIB_MGIS_CONFIG_HXX */
