/*!
 * \file   include/MFront/LibrariesManager.hxx
 * \brief
 * \author Thomas Helfer
 * \date   20/06/2018
 */

#ifndef LIB_MFRONT_LIBRARIESMANAGER_HXX
#define LIB_MFRONT_LIBRARIESMANAGER_HXX

#include <map>
#include <vector>
#include <string>

#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */

#include "MFront/Behaviour/Hypothesis.hxx"

namespace mfront {

  /*!
   * \brief structure in charge of handling libraries and querying MFront'
   * meta-data
   */
  struct LibrariesManager {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    //! a simple alias
    using libhandler = HINSTANCE__*;
#else
    //! a simple alias
    using libhandler = void*;
#endif /* LIB_EXTERNALLIBRARYMANAGER_HXX */
    //! a simple alias
    using Hypothesis = mfront::behaviour::Hypothesis;
    //! \return the only instance of this class (singleton pattern).
    static LibrariesManager& get();
    LibrariesManager(LibrariesManager&&) = delete;
    LibrariesManager(const LibrariesManager&) = delete;
    LibrariesManager& operator=(LibrariesManager&&) = delete;
    LibrariesManager& operator=(const LibrariesManager&) = delete;
    /*!
     * \return the `TFEL` version used to generate an `MFront` entry point
     * \param[in] l: library name
     * \param[in] n: entry point name
     */
    std::string getTFELVersion(const std::string&, const std::string&);
    /*!
     * \return the names of the material properties associated with a behaviour
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<std::string> getMaterialPropertiesNames(const std::string&,
                                                        const std::string&,
                                                        const Hypothesis);

   private:
    //! \brief constructor
    LibrariesManager();
    //! \brief destructor
    ~LibrariesManager();
    /*!
     * \return the names of a group of variables associated with a behaviour.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     * \param[in] n: name of the group of variables (such as MaterialProperties)
     */
    std::vector<std::string> getNames(const std::string&,
                                      const std::string&,
                                      const Hypothesis,
                                      const std::string&);
    /*!
     * \brief load the given symbol and cast it to a pointer of the given
     * type.
     * \tparam T: expected type
     * \param[in] l: library name
     * \param[in] n: symbol name
     */
    template <typename T>
    const T* extract(const std::string&, const std::string&);
    /*!
     * \brief load the given symbol and cast it to a pointer of the given
     * type.
     * \tparam T: expected type
     * \param[in] l: library name
     * \param[in] n1: symbol name
     * \param[in] n2: symbol name
     */
    template <typename T>
    const T* extract(const std::string&,
                     const std::string&,
                     const std::string&);
    /*!
     * \return the adress of the given symbol. If no symbol is found, a
     * null pointer is returned.
     * \param[in] l: library name
     * \param[in] n: symbol name
     */
    void* getSymbolAddress(const std::string&, const std::string&);
    /*!
     * \return the adress of the first given symbol if it exists or the
     * adress of the second symbol. If no symbol is found, a null pointer is
     * returned.
     * \param[in] l: library name
     * \param[in] n1: symbol name
     * \param[in] n2: symbol name
     */
    void* getSymbolAddress(const std::string&,
                           const std::string&,
                           const std::string&);
    /*!
     * \brief load an external library if not already loaded. If the library
     * is
     * successfully loaded, the associated handler is stored.
     * \param[in] l: library name
     * \return the handler to the library
     */
    libhandler loadLibrary(const std::string&);
    //! list of alreay loaded libraries
    std::map<std::string, libhandler> libraries;

    };  // end of struct LibrariesManager

}  // end of namespace mfront

#endif /* LIB_MFRONT_LIBRARIESMANAGER_HXX */
