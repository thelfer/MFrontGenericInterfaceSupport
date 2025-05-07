/*!
 * \file   MGIS/Context.hxx
 * \brief  This file declares the `Context` class
 * \date   08/02/2023
 */

#ifndef LIB_MGIS_CONTEXT_HXX
#define LIB_MGIS_CONTEXT_HXX 1

#include <ostream>
#include "MGIS/Config.hxx"
#include "MGIS/VerbosityLevel.hxx"
#include "MGIS/ErrorBacktrace.hxx"

namespace mgis {

  /*!
   * \brief a class used to pass an execution context to most methods of
   * `MGIS` and gather information (error, logs).
   *
   * The `Context` may be changed at various stage of the computation.
   * For example, the verbosity level or the logging stream
   * can be changed when calling a new model: this is can useful
   * to debug a specific rm.
   *
   * The default logging stream is the one returned by the
   * `mgis::getLogStream` free function.
   */
  class MGIS_EXPORT Context : public ErrorBacktrace {
   public:
    //! \brief default constructor
    Context() noexcept;
    //! \return the verbosity level
    [[nodiscard]] const VerbosityLevel &getVerbosityLevel() const noexcept;
    /*!
     * \brief change the level of verbosity
     * \param[in] l: the new verbose level
     */
    void setVerbosityLevel(const VerbosityLevel) noexcept;
    //! \brief destructor
    ~Context() noexcept;

   private:
    //
    Context(Context &&) = delete;
    Context(const Context &) = delete;
    Context &operator=(Context &&) = delete;
    Context &operator=(const Context &) = delete;
    /*!
     * \brief local level of verbosity, initialize by the global option
     * returned by the `getVerbosityLevel` function
     */
    VerbosityLevel verbosity_;
  };  // end of class Context

}  // end of namespace mgis

#endif /* LIB_MGIS_CONTEXT_HXX */
