/*!
 * \file   include/MGIS/Raise.ixx
 * \brief
 * \author Thomas Helfer
 * \date   20/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_RAISE_IXX
#define LIB_MGIS_RAISE_IXX

namespace mgis {

  template <typename Exception>
  void raise() {
    Exception e;
    throw(std::move(e));
  }  // end of raise

  template <typename Exception, typename... Args>
  void raise(Args&&... a) {
    Exception e(std::forward<Args...>(a...));
    throw(std::move(e));
  }  // end of raise

  template <typename Exception>
  void raise_if(const bool c) {
    if (c) {
      raise<Exception>();
    }
  }  // end of raise

  template <typename Exception, typename... Args>
  void raise_if(const bool c, Args&&... a) {
    if (c) {
      raise<Exception>(std::forward<Args...>(a...));
    }
  }  // end of raise

}  // end of namespace mgis

#endif /* LIB_MGIS_RAISE_IXX */
