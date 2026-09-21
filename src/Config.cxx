/*!
 * \file   Config.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/09/2026
 */

#include <cstdlib>
#include <iostream>
#include <exception>
#include "MGIS/Config.hxx"

namespace mgis::internal {

  static TerminateHandler &getTerminateHandler() noexcept {
    static TerminateHandler h = +[](std::string_view msg) {
      std::cerr << "mgis default terminate handler called\n"  //
                << msg << '\n';
      std::terminate();
    };
    return h;
  }  // end of getTerminateHandler

}  // end of namespace mgis::internal

namespace mgis {

  void setTerminateHandler(TerminateHandler &h) noexcept {
    if (h == nullptr) {
      return;
    }
    ::mgis::internal::getTerminateHandler() = h;
  }  // end of setTerminateHandler

  void terminate(std::string_view msg) {
    ::mgis::internal::getTerminateHandler()(msg);
    std::abort();  // just in case the handler returns
  }                // end of terminate

  void abort(std::string_view msg) {
    ::mgis::terminate(msg);
  }  // end of abort

}  // end of namespace mgis
