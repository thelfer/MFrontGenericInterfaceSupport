/*!
 * \file   Hypothesis.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_HYPOTHESIS_HXX
#define LIB_MGIS_BEHAVIOUR_HYPOTHESIS_HXX

#include <string>
#include "MGIS/Config.hxx"

namespace mgis {

  namespace behaviour {

    //! \brief the list of supported modelling hypotheses
    enum struct Hypothesis {
      AXISYMMETRICALGENERALISEDPLANESTRAIN,
      AXISYMMETRICALGENERALISEDPLANESTRESS,
      AXISYMMETRICAL,
      PLANESTRESS,
      PLANESTRAIN,
      GENERALISEDPLANESTRAIN,
      TRIDIMENSIONAL
    };  // end of enum Hypothesis

    /*!
     * \return the dimension of the space associated to the given hypothesis.
     * \param[in] h: modelling hypothesis
     */
    MGIS_EXPORT size_type getSpaceDimension(const Hypothesis);
    /*!
     * \return the size of a symmetric tensor in the given hypothesis.
     * \param[in] h: modelling hypothesis
     */
    MGIS_EXPORT size_type getStensorSize(const Hypothesis);
    /*!
     * \return the size of a tensor in the given hypothesis
     * \param[in] h: modelling hypothesis
     */
    MGIS_EXPORT size_type getTensorSize(const Hypothesis);
    /*!
     * \return the string associated to the given hypothesis
     * \param[in] h: modelling hypothesis
     */
    MGIS_EXPORT const char* toString(const Hypothesis);
    /*!
     * \return the hypothesis described by the given string.
     * \param[in] h: modelling hypothesis
     * \note valid values are:
     * - `AxisymmetricalGeneralisedPlaneStrain`
     * - `AxisymmetricalGeneralisedPlaneStress`
     * - `Axisymmetrical`
     * - `PlaneStress`
     * - `PlaneStrain`
     * - `GeneralisedPlaneStrain`
     * - `Tridimensional`
     */
    MGIS_EXPORT Hypothesis fromString(const char* const);
    /*!
     * \return the hypothesis described by the given string.
     * \param[in] h: modelling hypothesis
     * \note valid values are:
     * - `AxisymmetricalGeneralisedPlaneStrain`
     * - `AxisymmetricalGeneralisedPlaneStress`
     * - `Axisymmetrical`
     * - `PlaneStress`
     * - `PlaneStrain`
     * - `GeneralisedPlaneStrain`
     * - `Tridimensional`
     */
    MGIS_EXPORT Hypothesis fromString(const std::string&);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_HYPOTHESIS_HXX */
