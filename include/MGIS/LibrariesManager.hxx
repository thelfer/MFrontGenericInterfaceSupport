/*!
 * \file   include/MGIS/LibrariesManager.hxx
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

#ifndef LIB_MGIS_LIBRARIESMANAGER_HXX
#define LIB_MGIS_LIBRARIESMANAGER_HXX

#include <map>
#include <string>
#include <vector>

#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */

#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis {

/*!
 * \brief structure in charge of handling libraries and querying MGIS'
 * meta-data
 */
struct LibrariesManager {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
  //! a simple alias
  using libhandler = HINSTANCE__ *;
#else
  //! a simple alias
  using libhandler = void *;
#endif /* LIB_EXTERNALLIBRARYMANAGER_HXX */
  //! a simple alias
  using Hypothesis = mgis::behaviour::Hypothesis;
  //! \return the only instance of this class (singleton pattern).
  static LibrariesManager &get();
  LibrariesManager(LibrariesManager &&) = delete;
  LibrariesManager(const LibrariesManager &) = delete;
  LibrariesManager &operator=(LibrariesManager &&) = delete;
  LibrariesManager &operator=(const LibrariesManager &) = delete;
  /*!
   * \return the `TFEL` version used to generate an `MGIS` entry point
   * \param[in] l: library name
   * \param[in] n: entry point name
   */
  std::string getTFELVersion(const std::string &, const std::string &);
  /*!
   * \return the material knowledge type of an entry point in a
   * library. The returned value has the following meaning:
   * - 0: material property
   * - 1: behaviour
   * - 2: model
   * \param[in] l: library name
   */
  unsigned short getMaterialKnowledgeType(const std::string &,
                                          const std::string &);
  /*!
   * \return the name of the interface used to generate the given entry point
   * \param[in] l: library name
   * \param[in] n: entry point name
   */
  std::string getInterface(const std::string &, const std::string &);
  /*!
   * \return the name of the `MGIS` file used to generate the given entry point
   * \param[in] l: library name
   * \param[in] n: entry point name
   */
  std::string getSource(const std::string &, const std::string &);
  /*!
   * \return the type of the behaviour
   * \see MechanicalBehaviourBase::BehaviourType
   * The value returned has the following meaning:
   * - 0: general behaviour
   * - 1: strain based behaviour
   * - 2: standard finite strain behaviour
   * - 3: cohesive zone model
   * \param[in] l : name of the library
   * \param[in] b : behaviour name
   */
  unsigned short getBehaviourType(const std::string &, const std::string &);
  /*!
   * \return the kinematic assumption used by the behaviour
   * \see MechanicalBehaviourBase::Kinematic
   * The value returned has the following meaning:
   * - 0: undefined kinematic
   * - 1: standard small strain behaviour kinematic
   * - 2: cohesive zone model kinematic
   * - 3: standard finite strain kinematic (F-Cauchy)
   * - 4: ptest finite strain kinematic (eto-pk1)
   * - 5: Green-Lagrange strain
   * - 6: Miehe Apel Lambrecht logarithmic strain framework
   * \param[in] l : name of the library
   * \param[in] f : behaviour name
   */
  unsigned short getBehaviourKinematic(const std::string &,
                                       const std::string &);
  /*!
   * \return the symmetry of the behaviour
   * The value returned has the following meaning:
   * - 0: isotropic behaviour
   * - 1: orthotropic behaviour
   * \param[in] l : name of the library
   * \param[in] b : behaviour name
   */
  unsigned short getBehaviourSymmetry(const std::string &, const std::string &);
  /*!
   * \return the symmetry of the elastic stiffness
   * The value returned has the following meaning:
   * - 0: isotropic elastic stiffness
   * - 1: orthotropic elastic stiffness
   * \param[in] l : name of the library
   * \param[in] b : behaviour name
   */
  unsigned short getElasticStiffnessSymmetry(const std::string &,
                                             const std::string &);
  /*!
   * \return true if a behaviour generated throught the aster
   * interface requires a offset for the elastic properties
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  bool requiresStiffnessTensor(const std::string &, const std::string &,
                               const Hypothesis);
  /*!
   * \return true if a behaviour generated throught the aster
   * interface requires a offset for the elastic properties
   * \param[in] l: name of the library
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  bool requiresThermalExpansionCoefficientTensor(const std::string &,
                                                 const std::string &,
                                                 const Hypothesis);
  /*!
   * \return the names of the material properties associated with
   * a behaviour
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  std::vector<std::string> getMaterialPropertiesNames(const std::string &,
                                                      const std::string &,
                                                      const Hypothesis);
  /*!
   * \return the names of the internal state variables associated with
   * a behaviour
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  std::vector<std::string> getInternalStateVariablesNames(const std::string &,
                                                          const std::string &,
                                                          const Hypothesis);
  /*!
   * \return the types of the internal state variables associated with
   * a behaviour
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  std::vector<int> getInternalStateVariablesTypes(const std::string &,
                                                  const std::string &,
                                                  const Hypothesis);
  /*!
   * \return the names of the external state variables associated with
   * a behaviour
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   */
  std::vector<std::string> getExternalStateVariablesNames(const std::string &,
                                                          const std::string &,
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
   * \param[in] n: name of the group of variables (such as
   * MaterialProperties)
   */
  std::vector<std::string> getNames(const std::string &, const std::string &,
                                    const Hypothesis, const std::string &);
  /*!
   * \brief load the given symbol and cast it to a pointer of the given
   * type.
   * \tparam T: expected type
   * \param[in] l: library name
   * \param[in] n: symbol name
   */
  template <typename T>
  const T *extract(const std::string &, const std::string &);
  /*!
   * \brief load the given symbol and cast it to a pointer of the given
   * type.
   * \tparam T: expected type
   * \param[in] l: library name
   * \param[in] n1: symbol name
   * \param[in] n2: symbol name
   */
  template <typename T>
  const T *extract(const std::string &, const std::string &,
                   const std::string &);
  /*!
   * \return the adress of the given symbol. If no symbol is found, a
   * null pointer is returned.
   * \param[in] l: library name
   * \param[in] n: symbol name
   */
  void *getSymbolAddress(const std::string &, const std::string &);
  /*!
   * \return the adress of the first given symbol if it exists or the
   * adress of the second symbol. If no symbol is found, a null pointer is
   * returned.
   * \param[in] l: library name
   * \param[in] n1: symbol name
   * \param[in] n2: symbol name
   */
  void *getSymbolAddress(const std::string &, const std::string &,
                         const std::string &);
  /*!
   * \brief load an external library if not already loaded. If the library
   * is
   * successfully loaded, the associated handler is stored.
   * \param[in] l: library name
   * \return the handler to the library
   */
  libhandler loadLibrary(const std::string &);
  //! list of alreay loaded libraries
  std::map<std::string, libhandler> libraries;

}; // end of struct LibrariesManager

} // end of namespace mgis

#endif /* LIB_MGIS_LIBRARIESMANAGER_HXX */
