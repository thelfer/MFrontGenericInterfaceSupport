/*!
 * \file   MGIS/Context.hxx
 * \brief  This file declares the `Context` class
 * \date   08/02/2023
 */

#ifndef LIB_MGIS_CONTEXT_HXX
#define LIB_MGIS_CONTEXT_HXX 1

#include <memory>
#include <variant>
#include <ostream>
#include "MGIS/Config.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/LogStream.hxx"
#include "MGIS/VerbosityLevel.hxx"
#include "MGIS/ErrorBacktrace.hxx"

namespace mgis {

  /*!
   * \brief class that can be used to initialize a context, in particular in a
   * `constexpr` function.
   */
  struct ContextInitializer {
    //! \brief verbosity level
    const VerbosityLevel verbosity = VerbosityLevel::verboseQuiet;
  };

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
   * `mgis::getDefaultLogStream` free function.
   */
  struct MGIS_EXPORT Context final : public ErrorBacktrace {
    enum FailureHandlerPolicy { ABORT, RAISE };
    template <FailureHandlerPolicy policy>
    struct FailureHandler {
      FailureHandler(Context &c) : ctx(c) {}
      template <typename T>
      constexpr T &&operator()(std::optional<T> &&v) const {
        if (isInvalid(v)) {
          raise(this->ctx.getErrorMessage());
        }
        return std::move(*v);
      }
      template <typename T>
      constexpr T &operator()(OptionalReference<T> &&v) const {
        if (isInvalid(v)) {
          raise(this->ctx.getErrorMessage());
        }
        return *v;
      }

     private:
      //! \brief reference to the context that created the failure handler
      Context &ctx;
    };
    /*!
     * \brief default constructor
     *
     * The verbositiy level is initialized by calling the
     * `getDefaultVerbosityLevel` function.
     */
    Context() noexcept;
    /*!
     * \brief constructor for an initializer
     * \param[in] i: initializer
     */
    Context(const ContextInitializer &) noexcept;
    //! \return the verbosity level
    [[nodiscard]] const VerbosityLevel &getVerbosityLevel() const noexcept;
    /*!
     * \brief change the level of verbosity
     * \param[in] l: the new verbose level
     */
    void setVerbosityLevel(const VerbosityLevel) noexcept;
    /*!
     * \return a failure handler
     * \param[in] policy: policy used to treat a failure
     * \note the context must outlive the failure hander
     */
    template <FailureHandlerPolicy policy = FailureHandlerPolicy::RAISE>
    constexpr FailureHandler<policy> getFailureHandler() {
      return {*this};
    }
    /*!
     * \brief set the current log stream.
     * \param[in] s: log stream
     * \note the user is responsible for ensuring that the given object is alive
     */
    void setLogStream(std::ostream &) noexcept;
    /*!
     * \brief set the current log stream.
     * \param[in] s: log stream
     * \note if an empty shared pointer is given, the log stream is reset to the
     * default one, i.e. the one returned by the `mgis::getDefaultLogStream`
     * free function.
     */
    void setLogStream(std::shared_ptr<std::ostream>) noexcept;
    //! \brief reset the default log stream
    void resetLogStream() noexcept;
    /*!
     *  \brief disable the default log stream
     *
     * \note logging is disable by creating a no-op output stream
     */
    void disableLogStream() noexcept;
    /*!
     * \return the current log stream
     *
     * \note if no log stream is set, the default one is returned. See
     * `getDefaultLogStream` for details.
     */
    [[nodiscard]] std::ostream &log() noexcept;
    /*!
     * \brief display the given arguments in the log stream if the current
     * verbosity level (as returned by the `getVerbosityLevel` method) is
     * greater than a minimal one.
     *
     * \tparam Args: types of the arguments
     * \return the current log stream
     * \param[in] l: minimal verbosity level
     * \param[in] args: streamed object
     * \note note nothing is displayed if the current verbositiy
     * level is below the first argument
     */
    template <typename... Args>
    std::ostream &log(const VerbosityLevel, Args &&...) noexcept;
    /*!
     * \brief a simple wrapper around the `log` method to print a warning
     *
     * \tparam Args: types of the arguments
     * \param[in] args: streamed object
     */
    template <typename... Args>
    void warning(Args &&...) noexcept;
    /*!
     * \brief a simple wrapper around the `log` method which sets the minimun
     * verbosity level to `verboseDebug`
     *
     * \tparam Args: types of the arguments
     * \param[in] args: streamed object
     * \note note nothing is displayed if the current verbositiy level is below
     * `verboseDebug`
     */
    template <typename... Args>
    void debug(Args &&...) noexcept;
    //! \brief destructor
    ~Context() noexcept override;

   private:
    //
    Context(Context &&) = delete;
    Context(const Context &) = delete;
    Context &operator=(Context &&) = delete;
    Context &operator=(const Context &) = delete;
    //! \brief current log stream
    std::variant<std::monostate, std::ostream *, std::shared_ptr<std::ostream>>
        log_stream;
    /*!
     * \brief local level of verbosity, initialize by the
     * global option returned by the `getVerbosityLevel`
     * function
     */
    VerbosityLevel verbosity;
  };  // end of class Context

}  // end of namespace mgis

#include "MGIS/Context.ixx"

#endif /* LIB_MGIS_CONTEXT_HXX */
