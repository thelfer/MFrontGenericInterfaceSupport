/*!
 * \file   include/MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_FINITESTRAINBEHAVIOUROPTIONS_HXX
#define LIB_MGIS_BEHAVIOUR_FINITESTRAINBEHAVIOUROPTIONS_HXX

namespace mgis::behaviour {

  /*!
   * \brief option available for finite strain behaviours
   */
  struct FiniteStrainBehaviourOptions {
    //! \brief stress measure requested for finite strain behaviours
    enum StressMeasure {
      CAUCHY,  //!< Cauchy stress
      PK2,     //!< Second Piola-Kirchoff stress
      PK1      //!< First Piola-Kirchoff stress
    } stress_measure = CAUCHY;
    /*!
     * \brief type of finite strain tangent operator requested for finite
     * strain behaviours
     */
    enum TangentOperator {
      DSIG_DF, /*!< derivative of the Cauchy stress with respect to the
                    deformation gradient */
      DS_DEGL, /*!< derivative of the second Piola-Kirchoff stress with
                    respect to the Green-Lagrange strain */
      DPK1_DF, /*!< derivative of the first Piola-Kirchoff stress with
                    respect to the deformation gradient  */
      DTAU_DDF /*!< derivative of the Kirchoff stress with
                    respect to the spatial increment of the deformation gradient  */
    } tangent_operator = DSIG_DF;
  };  // end of struct FiniteStrainBehaviourOptions

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_FINITESTRAINBEHAVIOUROPTIONS_HXX */
