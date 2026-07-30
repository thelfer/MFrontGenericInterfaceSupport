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
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Raise.hxx"
#include "MGIS/LogStream.hxx"
#include "MGIS/VerbosityLevel.hxx"
#include "MGIS/ErrorBacktrace.hxx"
#include "MGIS/ProfilingData.hxx"
#include "MGIS/Profiling.hxx"

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
   * `MGIS` and gather information (error, logs, profiling).
   *
   * The default logging stream is the one returned by the
   * `mgis::getDefaultLogStream` free function.
   */
  struct MGIS_EXPORT Context final : public ErrorBacktrace {
    enum FailureHandlerPolicy { ABORT, RAISE };
    template <FailureHandlerPolicy policy>
    struct FailureHandler {
      explicit FailureHandler(Context &c) noexcept : ctx(c) {}
      template <typename T>
      decltype(auto) operator()(T &&v) const
          requires((internal::OptionalTraits<T>::isSpecialized) &&
                   (std::is_rvalue_reference_v<decltype(v)>)) {
        if (isInvalid(v)) {
          if (policy == FailureHandlerPolicy::RAISE) {
            raise(this->ctx.getErrorMessage());
          } else {
            this->ctx.abort();
          }
        }
        return internal::OptionalTraits<T>::deference(std::move(v));
      }
      template <typename T>
      decltype(auto) operator()(T &&v) const
          requires((!internal::OptionalTraits<T>::isSpecialized) &&
                   (std::is_rvalue_reference_v<decltype(v)>)) {
        if (isInvalid(v)) {
          raise(this->ctx.getErrorMessage());
        }
        return std::move(v);
      }
      template <typename T>
      friend decltype(auto) operator|(T &&v, const FailureHandler &h) requires(
          std::is_rvalue_reference_v<decltype(v)>) {
        return h(std::move(v));
      }  // end of operator|

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
    //
    Context(Context &&) = delete;
    Context(const Context &) = delete;
    Context &operator=(Context &&) = delete;
    Context &operator=(const Context &) = delete;
    //! \return the verbosity level
    [[nodiscard]] const VerbosityLevel &getVerbosityLevel() const noexcept;

    /*!
     * \brief change the level of verbosity
     * \param[in] l: the new verbose level
     */
    void setVerbosityLevel(const VerbosityLevel) noexcept;

    /*!
     * \brief enable or disable profiling
     * \param[in] b: a boolean value stating whether profiling shall be enabled
     *
     * \note when profiling is disabled, profiling sections introduce
     * almost no overhead.
     */
    void enableProfiling(const bool) noexcept;

    /*!
     * \return true if profiling is enabled, false otherwise
     */
    [[nodiscard]] bool isProfilingEnabled() const noexcept;

    /*!
     * \brief start a new profiling section
     * \param[in] name: name of the profiling section
     * \param[in] enabled: boolean stating whether this profiling section
     * shall effectively collect timing information
     *
     * \return a profiling section object
     *
     * \note the returned object relies on RAII semantics:
     * timing starts during construction and stops during destruction.
     */
    [[nodiscard]] ProfilingSection startNewProfiling(std::string,
                                                     bool) noexcept;

    /*!
     * \brief push a new profiling node into the current execution stack
     * \param[in] name: name of the profiling section
     *
     * \note this method is meant to be called internally by the
     * `Profiling` class during its construction.
     */
    void pushProfilingNode(std::string) noexcept;

    /*!
     * \brief pop the current profiling node from the execution stack and
     * accumulate time \param[in] dt: execution time of the section in seconds
     *
     * \note this method is meant to be called internally by the
     * `Profiling` class during its destruction.
     */
    void popProfilingNode(double) noexcept;

    /*!
     * \return the root node of the profiling results tree gathered during
     * execution
     */
    [[nodiscard]] const ProfilingData &getProfilingResultTree() const noexcept;

    /*!
     * \return a failure handler
     * \tparam policy: policy used to treat a failure
     * \note the context must outlive the failure hander
     */
    template <FailureHandlerPolicy policy = FailureHandlerPolicy::RAISE>
    [[nodiscard]] FailureHandler<policy> getFailureHandler() {
      return FailureHandler<policy>{*this};
    }

    //! \return a failure handler throwing exception in case of failure
    [[nodiscard]] FailureHandler<FailureHandlerPolicy::RAISE>
    getThrowingFailureHandler() noexcept;

    //! \return a failure handler aborting the execution in case of failure
    [[nodiscard]] FailureHandler<FailureHandlerPolicy::ABORT>
    getFatalFailureHandler() noexcept;

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

    //! \return a pointer to a log stream. This pointer may be null.
    [[nodiscard]] std::shared_ptr<std::ostream> getLogStreamPointer()
        const noexcept;

    //! \brief reset the default log stream
    void resetLogStream() noexcept;

    /*!
     * \brief disable the default log stream
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
    //! \brief printing the error message on the log stream and abort the
    //! execution
    [[noreturn]] void abort();

    //! \brief current log stream
    std::variant<std::monostate, std::ostream *, std::shared_ptr<std::ostream>>
        log_stream;

    /*!
     * \brief local level of verbosity, initialize by the
     * global option returned by the `getVerbosityLevel`
     * function
     */
    VerbosityLevel verbosity;

    //! \brief boolean stating whether profiling is enabled
    bool profiling_enabled = false;

    //! \brief root node of the profiling tree containing all recorded sections
    ProfilingData root_profiling_data;

    /*!
     * \brief Keeps track of the current path in the profiling tree.
     *
     * \note This stack does not manage memory. The lifetime and memory
     * management of the ProfilingData nodes are strictly handled by the tree
     * structure itself.
     */
    std::vector<ProfilingData *> profiling_stack;

  };  // end of class Context

}  // end of namespace mgis

#include "MGIS/Context.ixx"

#endif /* LIB_MGIS_CONTEXT_HXX */