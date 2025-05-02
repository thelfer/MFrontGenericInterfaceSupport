/*!
 * \file   MGIS/QuadratureFunction/MechanicalEvaluators.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_QUADRATUREFUNCTION_MECHANICALEVALUATORS_HXX
#define LIB_QUADRATUREFUNCTION_MECHANICALEVALUATORS_HXX

#ifdef MGIS_HAVE_TFEL
#include "MGIS/QuadratureFunction/MechanicalEvaluators/vonMisesStressEvaluator.hxx"
#include "MGIS/QuadratureFunction/MechanicalEvaluators/PrincipalStressEvaluator.hxx"
#include "MGIS/QuadratureFunction/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.hxx"
#endif MGIS_HAVE_TFEL

#endif /* LIB_QUADRATUREFUNCTION_MECHANICALEVALUATORS_HXX */
