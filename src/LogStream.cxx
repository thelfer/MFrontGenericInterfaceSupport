/*!
 * \file   LogStream.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#include <memory>
#include <variant>
#include <fstream>
#include <iostream>
#include "MGIS/LogStream.hxx"

namespace mgis {

  [[nodiscard]] static std::
      variant<std::monostate, std::ostream *, std::shared_ptr<std::ostream>>
          &getGlobalLogStream() noexcept {
    static std::variant<std::monostate, std::ostream *,
                        std::shared_ptr<std::ostream>>
        os;
    return os;
  }  // end of getGlobalLogStream

  std::ostream &getDefaultLogStream() noexcept {
    auto &os = getGlobalLogStream();
    if (std::holds_alternative<std::ostream *>(os)) {
      auto* const ptr = std::get<std::ostream *>(os);
      if (ptr == nullptr) {
        return std::cout;
      }
      return *ptr;
    } else if (std::holds_alternative<std::shared_ptr<std::ostream>>(os)) {
      auto ptr = std::get<std::shared_ptr<std::ostream>>(os);
      if (ptr.get() == nullptr) {
        return std::cout;
      }
      return *ptr;
    }
    return std::cout;
  }  // end of getDefaultLogStream

  std::pair<bool, std::string> setDefaultLogStream(
      std::string_view n) noexcept {
    auto ptr = std::make_shared<std::ofstream>(std::string{n});
    if (!(*ptr)) {
      return {false,
              "setDefaultLogStream: can't open file '" + std::string{n} + "'"};
    }
    setDefaultLogStream(ptr);
    return {true, ""};
  }  // end of setDefaultLogStream

  void setDefaultLogStream(std::ostream & os) noexcept{
    getGlobalLogStream() = &os;
  }  // end of setDefaultLogStream

  void setDefaultLogStream(std::shared_ptr<std::ostream> ptr) noexcept {
    if (ptr.get() == nullptr) {
      resetDefaultLogStream();
      return;
    }
    getGlobalLogStream() = ptr;
  }  // end of setDefaultLogStream

  void resetDefaultLogStream() noexcept {
    getGlobalLogStream() = std::monostate{};
  }  // end of resetDefaultLogStream

  void disableDefaultLogStream() noexcept {
    /*!
     * \brief a buffer which allows to create no-op output streams, i.e. streams
     * that does not nothing.
     *
     * see
     * https://stackoverflow.com/questions/11826554/standard-no-op-output-stream
     * for details
     */
    struct NullBuffer : public std::streambuf {
      int overflow(int c) override { return c; }
    };
    static NullBuffer null_buffer;
    static std::ostream null_stream(&null_buffer);
    setDefaultLogStream(null_stream);
  }  // end of disableDefaultLogStream


  /*!
   * \class  TerminalColors
   * \brief  This class contains char sequence corresponding to colors.
   * This enables to write ouput messages in color in the terminal.
   * \author Thomas Helfer
   * \date   26/07/2006
   */
  struct TerminalColors {
    /*!
     * \brief char sequence correponding to black.
     * \code
     * cout.write(TerminalColors::Black,sizeof(TerminalColors::Black));
     * \endcode
     */
    static constexpr char Black[5] = {033, '[', '3', '0', 'm'};
    /*!
     * \brief char sequence correponding to red.
     * \code
     * cout.write(TerminalColors::Red,sizeof(TerminalColors::Red));
     * \endcode
     */
    static constexpr char Red[5] = {033, '[', '3', '1', 'm'};
    /*!
     * \brief char sequence correponding to green.
     * \code
     * cout.write(TerminalColors::Green,sizeof(TerminalColors::Greeb));
     * \endcode
     */
    static constexpr char Green[5] = {033, '[', '3', '2', 'm'};
    /*!
     * \brief char sequence correponding to yellow.
     * \code
     * cout.write(TerminalColors::Yellow,sizeof(TerminalColors::Yellow));
     * \endcode
     */
    static constexpr char Yellow[5] = {033, '[', '3', '3', 'm'};
    /*!
     * \brief char sequence correponding to blue.
     * \code
     * cout.write(TerminalColors::Blue,sizeof(TerminalColors::Blue));
     * \endcode
     */
    static constexpr char Blue[5] = {033, '[', '3', '4', 'm'};
    /*!
     * \brief char sequence correponding to purple.
     * \code
     * cout.write(TerminalColors::Purple,sizeof(TerminalColors::Purple));
     * \endcode
     */
    static constexpr char Purple[5] = {033, '[', '3', '5', 'm'};
    /*!
     * \brief char sequence correponding to light blue.
     * \code
     * cout.write(TerminalColors::LightBlue,sizeof(TerminalColors::LightBlue));
     * \endcode
     */
    static constexpr char LightBlue[5] = {033, '[', '3', '6', 'm'};
    /*!
     * \brief char sequence correponding to white.
     * \code
     * cout.write(TerminalColors::White,sizeof(TerminalColors::White));
     * \endcode
     */
    static constexpr char White[5] = {033, '[', '3', '7', 'm'};
    /*!
     * \brief char sequence correponding to the reset command.
     * This causes the terminal to go back to its initial state.
     * \code
     * cout.write(TerminalColors::Reset,sizeof(TerminalColors::Reset));
     * \endcode
     */
    static constexpr char Reset[4] = {033, '[', 'm', 017};
  };

  void setErrorColor(std::ostream &os) noexcept {
    os.write(TerminalColors::Red, sizeof(TerminalColors::Red));
  }

  void setWarningColor(std::ostream &os) noexcept {
    os.write(TerminalColors::LightBlue, sizeof(TerminalColors::LightBlue));
  }

  void resetOutputColor(std::ostream & os) noexcept {
    os.write(TerminalColors::Reset, sizeof(TerminalColors::Reset));
  }

}  // end of namespace mgis