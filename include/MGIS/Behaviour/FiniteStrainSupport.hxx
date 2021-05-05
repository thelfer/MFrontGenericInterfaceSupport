/*!
 * \file   include/MGIS/Behaviour/FiniteStrainSupport.hxx
 * \brief
 * \author Thomas Helfer
 * \date   25/01/2019
 * \copyright (C) Copyright Thomas Helfer 2018-2019.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_FINITESTRAINSUPPORT_HXX
#define LIB_MGIS_BEHAVIOUR_FINITESTRAINSUPPORT_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Span.hxx"

namespace mgis {

  namespace behaviour {

    // forward declaration
    struct MaterialDataManager;
    // forward declaration
    struct BehaviourData;

    /*!
     * \brief list of the finite strain stress available
     */
    enum struct FiniteStrainStress {
      PK1  //! first Piola-Kirchhoff
    };     // end of enum struct FiniteStrainStress

    /*!
     * \brief list of the finite strain tangent operator available
     */
    enum struct FiniteStrainTangentOperator {
      DPK1_DF /*!< derivate of the first Piola-Kirchhoff with respect to the
               *   deformation gradient */
    };        // end of enum struct FiniteStrainTangentOperator

    /*!
     * \brief convert the Cauchy stress to the first Piola-Kirchhoff stress in
     * 2D. No bounds checks is made, use with care.
     * \note this does not work in plane stress, as the deformation gradient
     * does not provide the axial component.
     * \param[in] P: first Piola-Kirchhoff stress
     * \param[in] F: deformation gradient
     * \param[in] s: Cauchy stress
     */
    void convertFiniteStrainStress_PK1_2D(real* const,
                                          const real* const,
                                          const real* const);
    /*!
     * \brief convert the Cauchy stress to the first Piola-Kirchhoff stress in
     * 3D. No bounds checks is made, use with care.
     * \param[in] P: first Piola-Kirchhoff stress
     * \param[in] F: deformation gradient
     * \param[in] s: Cauchy stress
     */
    MGIS_EXPORT void convertFiniteStrainStress_PK1_3D(real* const,
                                                      const real* const,
                                                      const real* const);

    /*!
     * \param[out] s: new stress
     * \param[in] m: material data manager
     * \param[in] t: expected finite strain stress type
     */
    MGIS_EXPORT void convertFiniteStrainStress(mgis::span<real>&,
                                               const MaterialDataManager&,
                                               const FiniteStrainStress);
    /*!
     * \param[out] K: new tangent operator
     * \param[in] m: material data manager
     * \param[in] t: expected finite strain operator type
     */
    MGIS_EXPORT void convertFiniteStrainTangentOperator(
        mgis::span<real>&,
        const MaterialDataManager&,
        const FiniteStrainTangentOperator);
    /*!
     * \param[out] s: new stress
     * \param[in] d: behaviour data
     * \param[in] t: expected finite strain stress type
     */
    MGIS_EXPORT void convertFiniteStrainStress(mgis::span<real>&,
                                               const BehaviourData&,
                                               const FiniteStrainStress);
    /*!
     * \param[out] K: new tangent operator
     * \param[in] d: behaviour data
     * \param[in] t: expected finite strain operator type
     */
    MGIS_EXPORT void convertFiniteStrainTangentOperator(
        mgis::span<real>&,
        const BehaviourData&,
        const FiniteStrainTangentOperator);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_FINITESTRAINSUPPORT_HXX */
