/*!
 * \file   MGIS/LogStream.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#ifndef LIB_MGIS_LOGSTREAM_IXX
#define LIB_MGIS_LOGSTREAM_IXX

namespace mgis {

  template <typename... Args>
  void warning(std::ostream &os, Args &&...args) noexcept {
    os << "[warning]: ";
    (os << ... << std::forward<Args>(args));
    os << std::endl;
  }  // end of namespace mgis

  template <typename... Args>
  void debug(std::ostream &os, Args &&...args) noexcept {
    os << "[debug]: ";
    (os << ... << std::forward<Args>(args));
    os << std::endl;
  }  // end of namespace mgis

}  // end of namespace mgis

#endif /* LIB_MGIS_LOGSTREAM_IXX */
