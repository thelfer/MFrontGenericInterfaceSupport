/*!
 * \file   LogStream.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#ifndef LIB_MGIS_LOGSTREAM_HXX
#define LIB_MGIS_LOGSTREAM_HXX

#include <iosfwd>
#include <string>
#include <memory>
#include <utility>
#include <string_view>
#include "MGIS/Config.hxx"

namespace mgis {

  /*!
   * \return the default log stream
   *
   * \note by default, std::cout is returned
   */
  MGIS_EXPORT [[nodiscard]] std::ostream &getDefaultLogStream() noexcept;
  /*!
   * \brief set the default log stream as a file.
   *
   * \return a pair whose firt term indicates if the operation was successful.
   * In case of failure, the second term contains if an error message.
   *
   * \param[in] n: file name
   * \note In parallel, each process shall have its own output file or call
   * `disableDefaultLogStream` except on the root process
   */
  MGIS_EXPORT [[nodiscard]] std::pair<bool, std::string> setDefaultLogStream(
      std::string_view) noexcept;
  /*!
   * \brief set the default log stream from an exisiting output stream
   *
   * \param[in, out] os: output stream
   * \note the user is responsible for ensuring that the given object is alive
   */
  MGIS_EXPORT void setDefaultLogStream(std::ostream &) noexcept;
  /*!
   * \brief set the default log stream from an exisiting output stream
   *
   * \param[in, out] os: output stream
   * \note if the given pointer is null, `resetDefaultLogStream` is called
   */
  MGIS_EXPORT void setDefaultLogStream(std::shared_ptr<std::ostream>) noexcept;
  //! \brief reset the default log stream
  MGIS_EXPORT void resetDefaultLogStream() noexcept;
  /*!
   *  \brief disable the default log stream
   *
   * \note logging is disable by creating a no-op output stream
   */
  MGIS_EXPORT void disableDefaultLogStream() noexcept;
  /*!
   * \brief list of available colors for output stream
   * \note this only works on the terminal
   */
  enum struct OutputStreamColors {
    BLACK,
    RED,
    GREEN,
    YELLOW,
    BLUE,
    PURPLE,
    LIGHTBLUE,
    WHITE,
    RESET
  };
  /*!
   * \brief change the color of the output for warnings for red
   * \param[out] out: output stream
   * \param[in] c: color
   * \note this only works on the terminal
   */
  MGIS_EXPORT bool setStreamColor(std::ostream &,
                                  const OutputStreamColors) noexcept;
  /*!
   * \brief reset the color of the output
   * \param[out] out: output stream
   * \param[in] l: verbosity level
   */
  MGIS_EXPORT bool resetStreamColor(std::ostream &) noexcept;
  /*!
   * \brief print a warning message
   *
   * \tparam Args: types of the arguments
   * \param[in] args: streamed object
   */
  template <typename... Args>
  void warning(std::ostream &, Args &&...) noexcept;
  /*!
   * \brief print message for debugging
   *
   * \tparam Args: types of the arguments
   * \param[in] args: streamed object
   */
  template <typename... Args>
  void debug(std::ostream &, Args &&...) noexcept;

}  // end of namespace mgis

#include "MGIS/LogStream.ixx"

#endif /* LIB_MGIS_LOGSTREAM_HXX */
