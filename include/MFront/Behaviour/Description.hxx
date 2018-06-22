/*!
 * \file   Description.hxx
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

#ifndef LIB_MFRONT_BEHAVIOUR_DESCRIPTION_HXX
#define LIB_MFRONT_BEHAVIOUR_DESCRIPTION_HXX

#include <vector>
#include "MFront/Config.hxx"
#include "MFront/Behaviour/Variable.hxx"
#include "MFront/Behaviour/Behaviour.hxx"

namespace mfront {

  namespace behaviour {

    /*!
     * \brief structure describing a behaviour
     */
    struct MFRONT_EXPORT Description {
      //! \behaviour symmetry
      enum Symmetry { ISOTROPIC, ORTHOTROPIC };
      //! \brief constructor
      Description();
      //! \brief move constructor
      Description(Description &&);
      //! \brief copy constructor
      Description(const Description &);
      //! \brief move assignement
      Description &operator=(Description &&);
      //! \brief standard assignement
      Description &operator=(const Description &);
      //! \brief destructor
      ~Description();
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
      //! pointer to the function implementing the behaviour
      Behaviour b;
      //! \behaviour type
      enum {
        GENERALBEHAVIOUR,
        STANDARDSTRAINBASEDBEHAVIOUR,
        STANDARDFINITESTRAINBEHAVIOUR,
        COHESIVEZONEMODEL
      } btype;
      //! kinematic of the behaviour treated
      enum Kinematic {
        UNDEFINEDKINEMATIC,
        SMALLSTRAINKINEMATIC,
        COHESIVEZONEKINEMATIC,
        FINITESTRAINKINEMATIC_F_CAUCHY,
        FINITESTRAINKINEMATIC_ETO_PK1
      } kinematic;
      //! behaviour symmetry
      Symmetry symmetry;
      //! gradients
      std::vector<Variable> gradients;
      //! fluxes associated to gradients
      std::vector<Variable> fluxes;
      //! material properties
      std::vector<Variable> mps;
      //! internal state variables
      std::vector<Variable> isvs;
      //! external state variables
      std::vector<Variable> esvs;
      //! parameters
      std::vector<Variable> parameters;
    };  // end of struct Description

    /*!
     * \brief load the description of a behaviour from a library
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] h: modelling hypothesis
     * \return the behaviour description
     */
    MFRONT_EXPORT Description load(const std::string&,
                                   const std::string&,
                                   const Hypothesis);

  }  // end of namespace behaviour

}  // end of namespace mfront

#endif /* LIB_MFRONT_BEHAVIOUR_DESCRIPTION_HXX */

