/*!
 * \file   Markdown.cxx
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

#include <algorithm>
#include "MGIS/Utilities/Markdown.hxx"

namespace mgis::utilities {

  std::string get_heading_signs(const mgis::size_type l) {
    return std::string(std::max(l, mgis::size_type{1}), '#');
  }  // end of get_heading_signs

}  // end of namespace mgis::utilities
