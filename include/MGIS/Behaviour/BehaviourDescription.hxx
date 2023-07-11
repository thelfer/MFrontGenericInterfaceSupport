/*!
 * \file   include/MGIS/Behaviour/BehaviourDescription.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOURDESCRIPTION_HXX
#define LIB_MGIS_BEHAVIOURDESCRIPTION_HXX

#include <map>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Variable.hxx"

namespace mgis::behaviour {

  //! \brief structure describing a behaviour
  struct MGIS_EXPORT BehaviourDescription {
    /*!
     * \brief maximum number of behaviour options, whatever the kind of
     * behaviour treated.
     */
    static constexpr const mgis::size_type nopts = 2;
    //! \brief behaviour symmetry
    enum Symmetry { ISOTROPIC, ORTHOTROPIC };
    //! \brief name of the library providing the behaviour
    std::string library;
    //! \brief name of the behaviour
    std::string behaviour;
    //! \brief modelling hypothesis
    Hypothesis hypothesis;
    /*!
     * \brief name of the function (build using the behaviour name and the
     * modelling hypothesis)
     */
    std::string function;
    //! \brief name of the `MFront` source file used to generate the behaviour
    std::string source;
    //! \brief version of `TFEL` used to generate the behaviour
    std::string tfel_version;
    //! \brief unit system used by the behaviour
    std::string unit_system;
    //! \brief behaviour type
    enum BehaviourType {
      GENERALBEHAVIOUR,
      STANDARDSTRAINBASEDBEHAVIOUR,
      STANDARDFINITESTRAINBEHAVIOUR,
      COHESIVEZONEMODEL
    } btype;
    //! \brief kinematic of the behaviour treated
    enum Kinematic {
      UNDEFINEDKINEMATIC,
      SMALLSTRAINKINEMATIC,
      COHESIVEZONEKINEMATIC,
      FINITESTRAINKINEMATIC_F_CAUCHY,
      FINITESTRAINKINEMATIC_ETO_PK1
    } kinematic;
    //! \brief behaviour symmetry
    Symmetry symmetry;
    //! \brief gradients
    std::vector<Variable> gradients;
    //! \brief thermodynamic forces associated to gradients
    std::vector<Variable> thermodynamic_forces;
    //! \brief material properties
    std::vector<Variable> mps;
    //! \brief internal state variables
    std::vector<Variable> isvs;
    //! \brief external state variables
    std::vector<Variable> esvs;
    //! \brief tangent operator blocks
    std::vector<std::pair<Variable, Variable>> to_blocks;
    //! \brief real parameters
    std::vector<std::string> params;
    //! \brief integer parameters
    std::vector<std::string> iparams;
    //! \brief unsigned short parameters
    std::vector<std::string> usparams;
    /*!
     * \brief this boolean is true if the behaviour computes the
     * energy stored by the material per unit of volume in the reference
     * configuration. The physical meaning of this
     * energy depends on the behaviour considered.
     */
    bool computesStoredEnergy;
    /*!
     * \brief this boolean is true if the behaviour computes the
     * energy dissipated by the material per unit of volume in the
     * reference configuration. The physical meaning of this
     * energy depends on the behaviour considered.
     */
    bool computesDissipatedEnergy;
    //! \brief destructor
    ~BehaviourDescription();

   protected:
    //! \brief constructor
    BehaviourDescription();
    //! \brief move constructor
    BehaviourDescription(BehaviourDescription &&);
    //! \brief copy constructor
    BehaviourDescription(const BehaviourDescription &);
    //! \brief move assignement
    BehaviourDescription &operator=(BehaviourDescription &&);
    //! \brief standard assignement
    BehaviourDescription &operator=(const BehaviourDescription &);
  };  // end of struct BehaviourDescription

  /*!
   * \brief load the description of a behaviour from a library
   *
   * \param[out] d: behaviour description
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   * \return the behaviour description
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT void loadBehaviourDescription(BehaviourDescription &,
                                            const std::string &,
                                            const std::string &,
                                            const Hypothesis);

  /*!
   * \return the size of an array able to contain all the values of the
   * tangent operator
   * \param[in] b: behaviour
   */
  MGIS_EXPORT mgis::size_type getTangentOperatorArraySize(
      const BehaviourDescription &);
  /*!
   * \return if the given behaviour is a standard finite strain behaviour,
   * i.e. is a finite strain behaviour using the standard finite strain
   * kinematic
   * (called F-Cauchy although the stress measure can be chosen when
   * loading the behaviour)
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT bool isStandardFiniteStrainBehaviour(const std::string &,
                                                   const std::string &);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <typename T>
  T getParameterDefaultValue(const BehaviourDescription &, const std::string &);

  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT double getParameterDefaultValue<double>(
      const BehaviourDescription &, const std::string &);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT int getParameterDefaultValue<int>(const BehaviourDescription &,
                                                const std::string &);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT unsigned short getParameterDefaultValue<unsigned short>(
      const BehaviourDescription &, const std::string &);
  /*!
   * \return true if the given variable has bounds
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasBounds(const BehaviourDescription &, const std::string &);
  /*!
   * \return true if the given variable has a lower bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasLowerBound(const BehaviourDescription &,
                                 const std::string &);
  /*!
   * \return true if the given variable has a upper bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasUpperBound(const BehaviourDescription &,
                                 const std::string &);
  /*!
   * \return the lower bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getLowerBound(const BehaviourDescription &,
                                        const std::string &);
  /*!
   * \return the upper bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getUpperBound(const BehaviourDescription &,
                                        const std::string &);
  /*!
   * \return true if the given variable has bounds
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasPhysicalBounds(const BehaviourDescription &,
                                     const std::string &);
  /*!
   * \return true if the given variable has a lower physical bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasLowerPhysicalBound(const BehaviourDescription &,
                                         const std::string &);
  /*!
   * \return true if the given variable has a upper physical bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasUpperPhysicalBound(const BehaviourDescription &,
                                         const std::string &);
  /*!
   * \return the lower bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getLowerPhysicalBound(const BehaviourDescription &,
                                                const std::string &);
  /*!
   * \return the upper bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getUpperPhysicalBound(const BehaviourDescription &,
                                                const std::string &);
  /*!
   * \brief print a detailled (verbose) description of the behaviour using
   * a markdown format
   * \param[in] os: ouptut stream
   * \param[in] b: behaviour
   * \param[in] l: title level
   */
  MGIS_EXPORT void print_markdown(std::ostream &,
                                  const BehaviourDescription &,
                                  const mgis::size_type);
  /*!
   * \return the inputs of an initialize function
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] i: initialize function
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT std::vector<Variable> getBehaviourInitializeFunctionInputs(
      const std::string &,
      const std::string &,
      const std::string &,
      const Hypothesis);
  /*!
   * \return the outputs of a post-processing functions
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] i: post-processing
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT std::vector<Variable> getBehaviourPostProcessingOutputs(
      const std::string &,
      const std::string &,
      const std::string &,
      const Hypothesis);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOURDESCRIPTION_HXX */
