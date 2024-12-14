/*!
 * \file   src/Raise.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/04/2021
 */

#include "MGIS/Raise.hxx"

namespace mgis {

  static ExceptionHandler& getInternalExceptionHandler() {
    static ExceptionHandler h = nullptr;
    return h;
  }  // end of getInternalExceptionHandlerPointer

  void setExceptionHandler(ExceptionHandler h) {
    getInternalExceptionHandler() = h;
  }  // end of setExceptionHandler

  ExceptionHandler getExceptionHandler() {
    return getInternalExceptionHandler();
  }  // end of getExceptionHandler

}  // end of namespace mgis
