/*!
 * \file   MGIS/Function/MechanicalEvaluators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_FUNCTION_MECHANICALEVALUATORS_HXX
#define LIB_FUNCTION_MECHANICALEVALUATORS_HXX

#ifdef MGIS_HAVE_TFEL
#include "MGIS/Function/MechanicalEvaluators/vonMisesStressEvaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators/PrincipalStressEvaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators/HydrostaticStressEvaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.hxx"
#endif MGIS_HAVE_TFEL

#endif /* LIB_FUNCTION_MECHANICALEVALUATORS_HXX */
