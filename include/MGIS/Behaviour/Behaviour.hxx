/*!
 * \file   include/MGIS/Behaviour/Behaviour.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_DESCRIPTION_HXX
#define LIB_MGIS_BEHAVIOUR_DESCRIPTION_HXX

#include <iosfwd>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Span.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/RotationMatrix.hxx"
#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MGIS/Behaviour/BehaviourFctPtr.hxx"

namespace mgis::behaviour {

  /*!
   * \brief structure describing a behaviour
   */
  struct MGIS_EXPORT Behaviour {
    /*!
     * \brief maximum number of behaviour options, whatever the kind of
     * behaviour treated.
     */
    static constexpr const mgis::size_type nopts = 2;
    //! \brief behaviour symmetry
    enum Symmetry { ISOTROPIC, ORTHOTROPIC };
    //! \brief constructor
    Behaviour();
    //! \brief move constructor
    Behaviour(Behaviour &&);
    //! \brief copy constructor
    Behaviour(const Behaviour &);
    //! \brief move assignement
    Behaviour &operator=(Behaviour &&);
    //! \brief standard assignement
    Behaviour &operator=(const Behaviour &);
    //! \brief destructor
    ~Behaviour();
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
    //! \brief name of the `MGIS` source file used to generate the behaviour
    std::string source;
    //! \brief version of `TFEL` used to generate the behaviour
    std::string tfel_version;
    //! \brief pointer to the function implementing the behaviour
    BehaviourFctPtr b = nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of the gradients
     * from the global frame to the material frame.
     */
    RotateBehaviourGradientsFctPtr rotate_gradients_ptr = nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of an array of
     * gradients from the global frame to the material frame.
     */
    RotateArrayOfBehaviourGradientsFctPtr rotate_array_of_gradients_ptr =
        nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of the
     * thermodynamic forces from the material frame to the global frame.
     */
    RotateBehaviourThermodynamicForcesFctPtr rotate_thermodynamic_forces_ptr =
        nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of an array of
     * thermodynamic forces from the material frame to the global frame.
     */
    RotateArrayOfBehaviourThermodynamicForcesFctPtr
        rotate_array_of_thermodynamic_forces_ptr = nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of the tangent
     * operator blocks from the material frame to the global frame.
     */
    RotateBehaviourTangentOperatorBlocksFctPtr
        rotate_tangent_operator_blocks_ptr = nullptr;
    /*!
     * \brief pointer to a function implementing the rotation of an array of
     * tangent operator blocks from the material frame to the global frame.
     */
    RotateArrayOfBehaviourTangentOperatorBlocksFctPtr
        rotate_array_of_tangent_operator_blocks_ptr = nullptr;
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
    /*!
     * \brief behaviour options
     *
     * This is currently only meaningfull for finite strain
     * behaviours where the options stores:
     * - the stress measure used
     * - the consistent tangent operator expected
     *
     * For finite strain behaviours, options[1] holds the stress measure
     * used:
     * - if options[1] < 0.5, the Cauchy stress is used
     * - if 0.5 < options[1] < 1.5, the second Piola-Kirchoff stress is used
     * - if 1.5 < options[1] < 2.5, the first Piola-Kirchoff stress is used
     *
     * For finite strain behaviours, options[2] holds the consitent tangent
     * operator returned by the behaviour:
     * - if options[2]<0.5, the derivative of the Cauchy stress with respect
     *   to the deformation gradient is returned
     * - if 0.5<options[2]<1.5, the derivative of the second Piola-Kirchoff
     *   stress with respect to the Green-Lagrange strain
     *   is returned
     * - if 1.5<options[2]<2.5, the derivative of the first Piola-Kirchoff
     *   stress with respect to the deformation gradient is returned
     */
    std::vector<mgis::real> options;
  };  // end of struct Behaviour

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
   * \brief load the description of a behaviour from a library
   *
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   * \return the behaviour description
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT Behaviour load(const std::string &,
                             const std::string &,
                             const Hypothesis);
  /*!
   * \brief load the description of a finite strain behaviour from a library
   *
   * \note This method can also be used to load a finite strain behaviour.
   * In this case, the default options are used (the stress measure is Cauchy,
   * the tangent operator is the derivative of the Cauchy stress with respect
   * to the deformation gradient).
   *
   * \param[in] o: options
   * \param[in] l: library name
   * \param[in] b: behaviour name
   * \param[in] h: modelling hypothesis
   * \return the behaviour description
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT Behaviour load(const FiniteStrainBehaviourOptions &,
                             const std::string &,
                             const std::string &,
                             const Hypothesis);
  /*!
   * \return the size of an array able to contain all the values of the
   * tangent operator
   * \param[in] b: behaviour
   */
  MGIS_EXPORT mgis::size_type getTangentOperatorArraySize(const Behaviour &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out,in] g: gradients
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   * \note the rotation matrix argument can be given as a:
   * - an array of size 9 which is interpreted as a 3x3 rotation matrix.
   * - an array of size 9*n where n is the number of integration
   *   points which is interpreted as a field of 3x3 rotation matrix.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const mgis::span<const real> &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out] mg: array of gradients in the material frame
   * \param[in] b: behaviour description
   * \param[out] gg: array of gradients in the global frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   * \note the rotation matrix argument can be given as a:
   * - an array of size 9 which is interpreted as a 3x3 rotation matrix.
   * - an array of size 9*n where n is the number of integration
   *   points which is interpreted as a field of 3x3 rotation matrix.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const mgis::span<const real> &,
                                   const mgis::span<const real> &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out,in] g: gradients
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const RotationMatrix2D&);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out] mg: array of gradients in the material frame
   * \param[in] b: behaviour description
   * \param[out] gg: array of gradients in the global frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const mgis::span<const real> &,
                                   const RotationMatrix2D&);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out,in] g: gradients
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const RotationMatrix3D&);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out] mg: array of gradients in the material frame
   * \param[in] b: behaviour description
   * \param[out] gg: array of gradients in the global frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(mgis::span<real>,
                                   const Behaviour &,
                                   const mgis::span<const real> &,
                                   const RotationMatrix3D&);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out,in] tf: thermodynamics forces
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const mgis::span<const real> &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out,in] tf: thermodynamics forces
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const RotationMatrix2D &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out,in] tf: thermodynamics forces
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const RotationMatrix3D &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out] gtf: thermodynamics forces in the global frame
   * \param[in] b: behaviour description
   * \param[in] mtf: thermodynamics forces in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const mgis::span<const real> &,
                                             const mgis::span<const real> &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out] gtf: thermodynamics forces in the global frame
   * \param[in] b: behaviour description
   * \param[in] mtf: thermodynamics forces in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const mgis::span<const real> &,
                                             const RotationMatrix2D&);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out] gtf: thermodynamics forces in the global frame
   * \param[in] b: behaviour description
   * \param[in] mtf: thermodynamics forces in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(mgis::span<real>,
                                             const Behaviour &,
                                             const mgis::span<const real> &,
                                             const RotationMatrix3D &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out,in] K: tangent operator blocks
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const mgis::span<const real> &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out,in] K: tangent operator blocks
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const RotationMatrix2D &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out,in] K: tangent operator blocks
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const RotationMatrix3D &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out] gK: tangent operator blocks in the global frame
   * \param[in] b: behaviour description
   * \param[in] mK: tangent operator blocks in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const mgis::span<const real> &,
                                               const mgis::span<const real> &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out] gK: tangent operator blocks in the global frame
   * \param[in] b: behaviour description
   * \param[in] mK: tangent operator blocks in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const mgis::span<const real> &,
                                               const RotationMatrix2D &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out] gK: tangent operator blocks in the global frame
   * \param[in] b: behaviour description
   * \param[in] mK: tangent operator blocks in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(mgis::span<real>,
                                               const Behaviour &,
                                               const mgis::span<const real> &,
                                               const RotationMatrix3D &);
  /*!
   * \brief set the value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   * \param[in] v: parameter value
   */
  MGIS_EXPORT void setParameter(const Behaviour &,
                                const std::string &,
                                const double);
  /*!
   * \brief set the value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   * \param[in] v: parameter value
   */
  MGIS_EXPORT void setParameter(const Behaviour &,
                                const std::string &,
                                const int);
  /*!
   * \brief set the value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   * \param[in] v: parameter value
   */
  MGIS_EXPORT void setParameter(const Behaviour &,
                                const std::string &,
                                const unsigned short);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <typename T>
  T getParameterDefaultValue(const Behaviour &, const std::string &);

  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT double getParameterDefaultValue<double>(const Behaviour &,
                                                      const std::string &);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT int getParameterDefaultValue<int>(const Behaviour &,
                                                const std::string &);
  /*!
   * \return the default value of a parameter
   * \param[in] b: behaviour description
   * \param[in] n: parameter name
   */
  template <>
  MGIS_EXPORT unsigned short getParameterDefaultValue<unsigned short>(
      const Behaviour &, const std::string &);
  /*!
   * \return true if the given variable has bounds
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasBounds(const Behaviour &, const std::string &);
  /*!
   * \return true if the given variable has a lower bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasLowerBound(const Behaviour &, const std::string &);
  /*!
   * \return true if the given variable has a upper bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasUpperBound(const Behaviour &, const std::string &);
  /*!
   * \return the lower bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getLowerBound(const Behaviour &, const std::string &);
  /*!
   * \return the upper bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getUpperBound(const Behaviour &, const std::string &);
  /*!
   * \return true if the given variable has bounds
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasPhysicalBounds(const Behaviour &, const std::string &);
  /*!
   * \return true if the given variable has a lower physical bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasLowerPhysicalBound(const Behaviour &,
                                         const std::string &);
  /*!
   * \return true if the given variable has a upper physical bound
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT bool hasUpperPhysicalBound(const Behaviour &,
                                         const std::string &);
  /*!
   * \return the lower bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getLowerPhysicalBound(const Behaviour &,
                                                const std::string &);
  /*!
   * \return the upper bound of the given variable
   * \param[in] b: behaviour
   * \param[in] v: variable name
   */
  MGIS_EXPORT long double getUpperPhysicalBound(const Behaviour &,
                                                const std::string &);
  /*!
   * \brief print a detailled (verbose) description of the behaviour using
   * a markdown format
   * \param[in] os: ouptut stream
   * \param[in] b: behaviour
   * \param[in] l: title level
   */
  MGIS_EXPORT void print_markdown(std::ostream &,
                                  const Behaviour &,
                                  const mgis::size_type);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_DESCRIPTION_HXX */
