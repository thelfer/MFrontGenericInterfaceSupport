/*!
 * \file   MGIS/Databas.hxx
 * \brief
 * \author Thomas Helfer
 * \date   12/11/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_DATABASE_HXX
#define LIB_MGIS_DATABASE_HXX

#ifndef MGIS_HAVE_TFEL
#error "TFEL support is required to include this header"
#endif /* MGIS_HAVE_TFEL */

#include "MFront/MFrontDatabase.hxx"
#include "MGIS/Config.hxx"

namespace mgis {

  //! \return a global instance of a MFrontDatabase
  MGIS_VISIBILITY_EXPORT mfront::MFrontDatabase& getDatabase() noexcept;

}  // namespace mgis

#endif /* LIB_MGIS_DATABASE_HXX */
