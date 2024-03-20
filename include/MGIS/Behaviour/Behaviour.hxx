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

#ifndef LIB_MGIS_BEHAVIOUR_HXX
#define LIB_MGIS_BEHAVIOUR_HXX

#include <map>
#include <span>
#include <iosfwd>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/RotationMatrix.hxx"
#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MGIS/Behaviour/BehaviourFctPtr.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"

namespace mgis::behaviour {

  //! \brief structure describing an initialize function of a behaviour
  struct MGIS_EXPORT BehaviourInitializeFunction {
    //! \brief constructor
    BehaviourInitializeFunction();
    //! \brief move constructor
    BehaviourInitializeFunction(BehaviourInitializeFunction &&);
    //! \brief copy constructor
    BehaviourInitializeFunction(const BehaviourInitializeFunction &);
    //! \brief move assignement
    BehaviourInitializeFunction &operator=(BehaviourInitializeFunction &&);
    //! \brief standard assignement
    BehaviourInitializeFunction &operator=(const BehaviourInitializeFunction &);
    //! \brief destructor
    ~BehaviourInitializeFunction();
    //! \brief pointer to the initialize function
    BehaviourInitializeFctPtr f;
    //! \brief inputs of the initialize function
    std::vector<Variable> inputs;
  };  // end of BehaviourInitializeFunction

  //! \brief structure describing a post-processing of a behaviour
  struct MGIS_EXPORT BehaviourPostProcessing {
    //! \brief constructor
    BehaviourPostProcessing();
    //! \brief move constructor
    BehaviourPostProcessing(BehaviourPostProcessing &&);
    //! \brief copy constructor
    BehaviourPostProcessing(const BehaviourPostProcessing &);
    //! \brief move assignement
    BehaviourPostProcessing &operator=(BehaviourPostProcessing &&);
    //! \brief standard assignement
    BehaviourPostProcessing &operator=(const BehaviourPostProcessing &);
    //! \brief destructor
    ~BehaviourPostProcessing();
    //! \brief pointer to the post-processing
    BehaviourPostProcessingFctPtr f;
    //! \brief outputs of the post-processing
    std::vector<Variable> outputs;
  };  // end of BehaviourPostProcessing

  //! \brief structure describing a behaviour
  struct MGIS_EXPORT Behaviour : BehaviourDescription {
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
    //! \brief list of initialize functions associated with the behaviour
    std::map<std::string, BehaviourInitializeFunction, std::less<>>
        initialize_functions;
    //! \brief pointer to the function implementing the behaviour
    BehaviourFctPtr b = nullptr;
    //! \brief list of post-processings associated with the behaviour
    std::map<std::string, BehaviourPostProcessing, std::less<>> postprocessings;
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
    /*!
     * \brief behaviour options
     *
     * This is currently only meaningfull for finite strain
     * behaviours where the options stores:
     * - the stress measure used
     * - the consistent tangent operator expected
     *
     * For finite strain behaviours, options[0] holds the stress measure
     * used:
     * - if options[0] < 0.5, the Cauchy stress is used
     * - if 0.5 < options[0] < 1.5, the second Piola-Kirchoff stress is used
     * - if 1.5 < options[0] < 2.5, the first Piola-Kirchoff stress is used
     *
     * For finite strain behaviours, options[1] holds the consitent tangent
     * operator returned by the behaviour:
     * - if options[1]<0.5, the derivative of the Cauchy stress with respect
     *   to the deformation gradient is returned
     * - if 0.5<options[1]<1.5, the derivative of the second Piola-Kirchoff
     *   stress with respect to the Green-Lagrange strain
     *   is returned
     * - if 1.5<options[1]<2.5, the derivative of the first Piola-Kirchoff
     *   stress with respect to the deformation gradient is returned
     */
    std::vector<mgis::real> options;
  };  // end of struct Behaviour

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
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const std::span<const real> &);
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
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const std::span<const real> &,
                                   const std::span<const real> &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out,in] g: gradients
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const RotationMatrix2D &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out] mg: array of gradients in the material frame
   * \param[in] b: behaviour description
   * \param[out] gg: array of gradients in the global frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const std::span<const real> &,
                                   const RotationMatrix2D &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out,in] g: gradients
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const RotationMatrix3D &);
  /*!
   * \brief rotate an array of gradients from the global frame to the material
   * frame.
   * \param[out] mg: array of gradients in the material frame
   * \param[in] b: behaviour description
   * \param[out] gg: array of gradients in the global frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateGradients(std::span<real>,
                                   const Behaviour &,
                                   const std::span<const real> &,
                                   const RotationMatrix3D &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out,in] tf: thermodynamics forces
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
                                             const Behaviour &,
                                             const std::span<const real> &);
  /*!
   * \brief rotate an array of thermodynamics forces from the material frame
   * to
   * the global frame.
   * \param[out,in] tf: thermodynamics forces
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
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
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
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
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
                                             const Behaviour &,
                                             const std::span<const real> &,
                                             const std::span<const real> &);
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
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
                                             const Behaviour &,
                                             const std::span<const real> &,
                                             const RotationMatrix2D &);
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
  MGIS_EXPORT void rotateThermodynamicForces(std::span<real>,
                                             const Behaviour &,
                                             const std::span<const real> &,
                                             const RotationMatrix3D &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out,in] K: tangent operator blocks
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
                                               const Behaviour &,
                                               const std::span<const real> &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out,in] K: tangent operator blocks
   * \param[in] b: behaviour description
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
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
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
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
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
                                               const Behaviour &,
                                               const std::span<const real> &,
                                               const std::span<const real> &);
  /*!
   * \brief rotate an array of tangent operator blocks from the material frame
   * to the global frame.
   * \param[out] gK: tangent operator blocks in the global frame
   * \param[in] b: behaviour description
   * \param[in] mK: tangent operator blocks in the material frame
   * \param[in] r: rotation matrix from the global frame to the material
   * frame.
   */
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
                                               const Behaviour &,
                                               const std::span<const real> &,
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
  MGIS_EXPORT void rotateTangentOperatorBlocks(std::span<real>,
                                               const Behaviour &,
                                               const std::span<const real> &,
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
   * \return the size of an array able to contain the inputs of an
   * initialize function.
   * \param[in] b: behaviour
   * \param[in] n: name of the post-processing
   */
  MGIS_EXPORT size_type getInitializeFunctionVariablesArraySize(
      const Behaviour &, const std::string_view);
  /*!
   * \return an array containing the inputs of an initialize function.
   * \param[in] b: behaviour
   * \param[in] n: name of the initialize function
   */
  MGIS_EXPORT std::vector<mgis::real> allocateInitializeFunctionVariables(
      const Behaviour &, const std::string_view);
  /*!
   * \return the size of an array able to contain the outputs of a
   * post-processing.
   * \param[in] b: behaviour
   * \param[in] n: name of the post-processing
   */
  MGIS_EXPORT size_type getPostProcessingVariablesArraySize(
      const Behaviour &, const std::string_view);
  /*!
   * \return an array containing the results of a post-processing.
   * \param[in] b: behaviour
   * \param[in] n: name of the post-processing
   */
  MGIS_EXPORT std::vector<mgis::real> allocatePostProcessingVariables(
      const Behaviour &, const std::string_view);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_HXX */
