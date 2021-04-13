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

#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MGIS/Behaviour/BehaviourFctPtr.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis {

  /*!
   * \brief structure in charge of handling libraries and querying MGIS'
   * meta-data
   */
  struct MGIS_VISIBILITY_EXPORT LibrariesManager {
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
     * \return the name of the `MGIS` file used to generate the given entry
     * point
     * \param[in] l: library name
     * \param[in] n: entry point name
     */
    std::string getSource(const std::string &, const std::string &);
    /*!
     * \return the function implementing the behaviour
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::BehaviourFctPtr getBehaviour(const std::string &,
                                                  const std::string &,
                                                  const Hypothesis);
    /*!
     * \return the function implementing the rotation of the gradients of a
     * behaviour from the global frame to the material frame
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateBehaviourGradientsFctPtr
    getRotateBehaviourGradientsFunction(const std::string &,
                                        const std::string &,
                                        const Hypothesis);
    /*!
     * \return the function implementing the rotation of the gradients of a
     * behaviour from the global frame to the material frame
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateArrayOfBehaviourGradientsFctPtr
    getRotateArrayOfBehaviourGradientsFunction(const std::string &,
                                               const std::string &,
                                               const Hypothesis);
    /*!
     * \return the function implementing the rotation of the thermodynamic
     * forces from the material frame to the global frame.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr
    getRotateBehaviourThermodynamicForcesFunction(const std::string &,
                                                  const std::string &,
                                                  const Hypothesis);
    /*!
     * \return the function implementing the rotation of the thermodynamic
     * forces from the material frame to the global frame.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr
    getRotateArrayOfBehaviourThermodynamicForcesFunction(const std::string &,
                                                         const std::string &,
                                                         const Hypothesis);
    /*!
     * \return the function implementing the rotation of the thermodynamic
     * forces from the material frame to the global frame for a finite strain
     * behaviour.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     * \param[in] s: stress returned by the behaviour
     */
    mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr
    getRotateBehaviourThermodynamicForcesFunction(
        const std::string &,
        const std::string &,
        const Hypothesis,
        const mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure);
    /*!
     * \return the function implementing the rotation of the thermodynamic
     * forces from the material frame to the global frame for a finite strain
     * behaviour.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     * \param[in] s: stress returned by the behaviour
     */
    mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr
    getRotateArrayOfBehaviourThermodynamicForcesFunction(
        const std::string &,
        const std::string &,
        const Hypothesis,
        const mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure);
    /*!
     * \return the function implementing the rotation of the tangent operator
     * block from the material frame to the global frame.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr
    getRotateBehaviourTangentOperatorBlocksFunction(const std::string &,
                                                    const std::string &,
                                                    const Hypothesis);
    /*!
     * \return the function implementing the rotation of the tangent operator
     * block from the material frame to the global frame.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     */
    mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr
    getRotateArrayOfBehaviourTangentOperatorBlocksFunction(const std::string &,
                                                           const std::string &,
                                                           const Hypothesis);
    /*!
     * \return the function implementing the rotation of the tangent operator
     * block from the material frame to the global frame for a finite strain
     * behaviour.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     * \param[in] t: tangent operator returned by the behaviour
     */
    mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr
    getRotateBehaviourTangentOperatorBlocksFunction(
        const std::string &,
        const std::string &,
        const Hypothesis,
        const mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator);
    /*!
     * \return the function implementing the rotation of the tangent operator
     * block from the material frame to the global frame for a finite strain
     * behaviour.
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: hypothesis
     * \param[in] t: tangent operator returned by the behaviour
     */
    mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr
    getRotateArrayOfBehaviourTangentOperatorBlocksFunction(
        const std::string &,
        const std::string &,
        const Hypothesis,
        const mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator);
    /*!
     * \return the type of the behaviour
     * \see MechanicalBehaviourBase::BehaviourType
     * The value returned has the following meaning:
     * - 0: general behaviour
     * - 1: strain based behaviour
     * - 2: standard finite strain behaviour
     * - 3: cohesive zone model
     * \param[in] l: library
     * \param[in] b: behaviour name
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
     * \param[in] l: library
     * \param[in] f: behaviour name
     */
    unsigned short getBehaviourKinematic(const std::string &,
                                         const std::string &);
    /*!
     * \return the symmetry of the behaviour
     * The value returned has the following meaning:
     * - 0: isotropic behaviour
     * - 1: orthotropic behaviour
     * \param[in] l: library
     * \param[in] b: behaviour name
     */
    unsigned short getBehaviourSymmetry(const std::string &,
                                        const std::string &);
    /*!
     * \return the symmetry of the elastic stiffness
     * The value returned has the following meaning:
     * - 0: isotropic elastic stiffness
     * - 1: orthotropic elastic stiffness
     * \param[in] l: library
     * \param[in] b: behaviour name
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
    bool requiresStiffnessTensor(const std::string &,
                                 const std::string &,
                                 const Hypothesis);
    /*!
     * \return true if a behaviour computes the internal energy. In `MFront`,
     * this is associated with the `@InternalEnergy` code block.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    bool computesStoredEnergy(const std::string &,
                              const std::string &,
                              const Hypothesis);
    /*!
     * \return true if a behaviour computes the dissipated energy. In `MFront`,
     * this is associated with the `@DissipatedEnergy` code block.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    bool computesDissipatedEnergy(const std::string &,
                                  const std::string &,
                                  const Hypothesis);
    /*!
     * \return true if a behaviour generated throught the aster
     * interface requires a offset for the elastic properties
     * \param[in] l: library
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    bool requiresThermalExpansionCoefficientTensor(const std::string &,
                                                   const std::string &,
                                                   const Hypothesis);
    /*!
     * \return the names of the gradients associated with a behaviour
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<std::string> getGradientsNames(const std::string &,
                                               const std::string &,
                                               const Hypothesis);
    /*!
     * \return the types of the gradients associated with a behaviour
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<int> getGradientsTypes(const std::string &,
                                       const std::string &,
                                       const Hypothesis);
    /*!
     * \return the names of the thermodynamic forces associated with a
     * behaviour.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<std::string> getThermodynamicForcesNames(const std::string &,
                                                         const std::string &,
                                                         const Hypothesis);
    /*!
     * \return the names of the tangent operator blocks names.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<std::pair<std::string, std::string>>
    getTangentOperatorBlocksNames(const std::string &,
                                  const std::string &,
                                  const Hypothesis);
    /*!
     * \return the types of the thermodynamic forces associated with a
     * behaviour.
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     */
    std::vector<int> getThermodynamicForcesTypes(const std::string &,
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
    /*!
     * \param[in] l: library
     * \param[in] s: name of function or mechanical behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     * \param[in] v: value
     */
    void setParameter(const std::string &,
                      const std::string &,
                      const Hypothesis,
                      const std::string &,
                      const double);
    /*!
     * set the value of an integer parameter
     * \param[in] l: library
     * \param[in] s: name of function or mechanical behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     * \param[in] v: value
     */
    void setParameter(const std::string &,
                      const std::string &,
                      const Hypothesis,
                      const std::string &,
                      const int);
    /*!
     * set the value of an unsigned short parameter
     * \param[in] l: library
     * \param[in] s: name of function or mechanical behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     * \param[in] v: value
     */
    void setParameter(const std::string &,
                      const std::string &,
                      const Hypothesis,
                      const std::string &,
                      const unsigned short);
    /*!
     * \param[in] l: library
     * \param[in] f: law name
     * \param[in] h: modelling hypothesis
     */
    std::vector<std::string> getParametersNames(const std::string &,
                                                const std::string &,
                                                const Hypothesis);
    /*!
     * \return the types associated with each parameter. The integer
     * values returned have the following meaning:
     *
     * - 0: floatting point value
     * - 1: integer value
     * - 2: unsigned short value
     *
     * \param[in] l: library
     * \param[in] f: law name
     * \param[in] h: modelling hypothesis
     */
    std::vector<int> getParametersTypes(const std::string &,
                                        const std::string &,
                                        const Hypothesis);
    /*!
     * \brief get the default value of a double parameter
     * \param[in] l: library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     */
    double getParameterDefaultValue(const std::string &,
                                    const std::string &,
                                    const Hypothesis,
                                    const std::string &);
    /*!
     * \brief get the default value of an integer parameter
     * \param[in] l: library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     */
    int getIntegerParameterDefaultValue(const std::string &,
                                        const std::string &,
                                        const Hypothesis,
                                        const std::string &);
    /*!
     * \brief get the default value of an unsigned short parameter
     * \param[in] l: library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameter name
     */
    unsigned short getUnsignedShortParameterDefaultValue(const std::string &,
                                                         const std::string &,
                                                         const Hypothesis,
                                                         const std::string &);
    /*!
     * \return true if the given variable has bounds
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasBounds(const std::string &,
                   const std::string &,
                   const Hypothesis,
                   const std::string &);
    /*!
     * \return true if the given variable has a lower bound
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasLowerBound(const std::string &,
                       const std::string &,
                       const Hypothesis,
                       const std::string &);
    /*!
     * \return true if the given variable has a upper bound
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasUpperBound(const std::string &,
                       const std::string &,
                       const Hypothesis,
                       const std::string &);
    /*!
     * \return the lower bound of the given variable
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    long double getLowerBound(const std::string &,
                              const std::string &,
                              const Hypothesis,
                              const std::string &);
    /*!
     * \return the upper bound of the given variable
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    long double getUpperBound(const std::string &,
                              const std::string &,
                              const Hypothesis,
                              const std::string &);
    /*!
     * \return true if the given variable has bounds
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasPhysicalBounds(const std::string &,
                           const std::string &,
                           const Hypothesis,
                           const std::string &);
    /*!
     * \return true if the given variable has a lower physical bound
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasLowerPhysicalBound(const std::string &,
                               const std::string &,
                               const Hypothesis,
                               const std::string &);
    /*!
     * \return true if the given variable has a upper physical bound
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    bool hasUpperPhysicalBound(const std::string &,
                               const std::string &,
                               const Hypothesis,
                               const std::string &);
    /*!
     * \return the lower bound of the given variable
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    long double getLowerPhysicalBound(const std::string &,
                                      const std::string &,
                                      const Hypothesis,
                                      const std::string &);
    /*!
     * \return the upper bound of the given variable
     * \param[in] l: name of the library
     * \param[in] b: behaviour
     * \param[in] h: modelling hypothesis
     * \param[in] v: variable name
     */
    long double getUpperPhysicalBound(const std::string &,
                                      const std::string &,
                                      const Hypothesis,
                                      const std::string &);

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
    std::vector<std::string> getNames(const std::string &,
                                      const std::string &,
                                      const Hypothesis,
                                      const std::string &);
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
    const T *extract(const std::string &,
                     const std::string &,
                     const std::string &);
    /*!
     * \return the given symbol exists.
     * \param[in] l: library name
     * \param[in] n: symbol name
     */
    bool contains(const std::string &, const std::string &);
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
    void *getSymbolAddress(const std::string &,
                           const std::string &,
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

  };  // end of struct LibrariesManager

}  // end of namespace mgis

#endif /* LIB_MGIS_LIBRARIESMANAGER_HXX */
