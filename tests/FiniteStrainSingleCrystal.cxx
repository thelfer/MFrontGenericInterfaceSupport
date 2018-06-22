/*!
 * \file   FiniteStrainSingleCrystal.cxx
 * \brief
 * \author Thomas Helfer
 * \date   22/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#define MFRONT_SHAREDOBJ __declspec(dllexport)
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (!defined __INTEL_COMPILER)
#if __GNUC__ >= 4
#define MFRONT_SHAREDOBJ __attribute__((visibility("default")))
#else /*__GNUC__ >= 4 */
#define MFRONT_SHAREDOBJ
#endif /*__GNUC__ >= 4 */
#elif defined __INTEL_COMPILER
#define MFRONT_SHAREDOBJ __attribute__((visibility("default")))
#else /* (defined __GNUC__) && (! defined __INTEL_COMPILER) */
#define MFRONT_SHAREDOBJ
#endif /* (defined __GNUC__) && (! defined __INTEL_COMPILER) */
#endif /* defined _WIN32 || defined _WIN64 ||defined __CYGWIN__ */

extern "C" {

MFRONT_SHAREDOBJ const char* FiniteStrainSingleCrystal_mfront_ept =
    "FiniteStrainSingleCrystal";

MFRONT_SHAREDOBJ const char* FiniteStrainSingleCrystal_tfel_version =
    "3.2.0-dev";

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char* FiniteStrainSingleCrystal_mfront_interface =
    "generic";

MFRONT_SHAREDOBJ const char* FiniteStrainSingleCrystal_src =
    "FiniteStrainSingleCrystal.mfront";

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_nModellingHypotheses =
    1u;

MFRONT_SHAREDOBJ const char* FiniteStrainSingleCrystal_ModellingHypotheses[1u] =
    {"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_BehaviourType = 2u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_BehaviourKinematic =
    3u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_SymmetryType = 1u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_ElasticSymmetryType =
    1u;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_OrthotropyManagementPolicy = 2u;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_UsableInPurelyImplicitResolution =
        1;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nMaterialProperties = 7u;

MFRONT_SHAREDOBJ const char*
    FiniteStrainSingleCrystal_Tridimensional_MaterialProperties[7u] = {
        "m", "K", "C", "R0", "Q", "b", "d1"};

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nInternalStateVariables = 37;
MFRONT_SHAREDOBJ const char*
    FiniteStrainSingleCrystal_Tridimensional_InternalStateVariables[37] = {
        "g[0]",  "g[1]", "g[2]",  "g[3]",  "g[4]", "g[5]", "g[6]", "g[7]",
        "g[8]",  "g[9]", "g[10]", "g[11]", "Fe",   "p[0]", "p[1]", "p[2]",
        "p[3]",  "p[4]", "p[5]",  "p[6]",  "p[7]", "p[8]", "p[9]", "p[10]",
        "p[11]", "a[0]", "a[1]",  "a[2]",  "a[3]", "a[4]", "a[5]", "a[6]",
        "a[7]",  "a[8]", "a[9]",  "a[10]", "a[11]"};
MFRONT_SHAREDOBJ int
    FiniteStrainSingleCrystal_Tridimensional_InternalStateVariablesTypes[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char* const*
    FiniteStrainSingleCrystal_Tridimensional_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nParameters = 12;
MFRONT_SHAREDOBJ const char*
    FiniteStrainSingleCrystal_Tridimensional_Parameters[12] = {
        "theta",
        "epsilon",
        "h1",
        "h2",
        "h3",
        "h4",
        "h5",
        "h6",
        "minimal_time_step_scaling_factor",
        "maximal_time_step_scaling_factor",
        "numerical_jacobian_epsilon",
        "iterMax"};
MFRONT_SHAREDOBJ int
    FiniteStrainSingleCrystal_Tridimensional_ParametersTypes[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2};

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_theta_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_epsilon_ParameterDefaultValue =
        1e-11;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h1_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h2_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h3_ParameterDefaultValue = 0.6;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h4_ParameterDefaultValue = 12.3;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h5_ParameterDefaultValue = 1.6;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h6_ParameterDefaultValue = 1.8;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_minimal_time_step_scaling_factor_ParameterDefaultValue =
        0.1;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_maximal_time_step_scaling_factor_ParameterDefaultValue =
        1.79769e+308;

MFRONT_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_numerical_jacobian_epsilon_ParameterDefaultValue =
        1e-12;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_iterMax_ParameterDefaultValue =
        100;

MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_requiresStiffnessTensor = 1;
MFRONT_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_requiresThermalExpansionCoefficientTensor =
        0;

}  // end of extern "C"
