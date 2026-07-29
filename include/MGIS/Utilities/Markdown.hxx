/*!
 * \file   Markdown.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_UTILITIES_MARKDOWN_HXX
#define LIB_MGIS_UTILITIES_MARKDOWN_HXX

#include <string>
#include "MGIS/Config.hxx"

namespace mgis::utilities {

  /*!
   * \return the heading signs associated with the given heading level
   * \param[in] l: heading level
   */
  MGIS_EXPORT std::string get_heading_signs(const mgis::size_type);

}  // end of namespace mgis::utilities

#endif /* LIB_MGIS_UTILITIES_MARKDOWN_HXX */
