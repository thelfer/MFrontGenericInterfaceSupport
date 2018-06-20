/*!
 * \file   Description.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
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
    struct Description {
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
      } ktype;
      //! \brief constructor
      Description();
      //! \brief move constructor
      Description(Description&&);
      //! \brief copy constructor
      Description(const Description&);
      //! \brief move assignement
      Description& operator=(Description&&);
      //! \brief standard assignement
      Description& operator=(const Description&);
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
      //! pointer to the function implementing the behaviour
      Behaviour b;
      //! gradients
      std::vector<Variable> gradients;
      //! fluxes associated to gradients
      std::vector<Variable> fluxes;
      //! material properties
      std::vector<Variable> mps;
      //! internal state variables
      std::vector<Variable> ivs;
      //! external state variables
      std::vector<Variable> evs;
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
