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
#define MGIS_SHAREDOBJ __declspec(dllexport)
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (!defined __INTEL_COMPILER)
#if __GNUC__ >= 4
#define MGIS_SHAREDOBJ __attribute__((visibility("default")))
#else /*__GNUC__ >= 4 */
#define MGIS_SHAREDOBJ
#endif /*__GNUC__ >= 4 */
#elif defined __INTEL_COMPILER
#define MGIS_SHAREDOBJ __attribute__((visibility("default")))
#else /* (defined __GNUC__) && (! defined __INTEL_COMPILER) */
#define MGIS_SHAREDOBJ
#endif /* (defined __GNUC__) && (! defined __INTEL_COMPILER) */
#endif /* defined _WIN32 || defined _WIN64 ||defined __CYGWIN__ */

extern "C" {

MGIS_SHAREDOBJ const char* FiniteStrainSingleCrystal_mfront_ept =
    "FiniteStrainSingleCrystal";

MGIS_SHAREDOBJ const char* FiniteStrainSingleCrystal_tfel_version =
    "3.2.0-dev";

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_mfront_mkt = 1u;

MGIS_SHAREDOBJ const char* FiniteStrainSingleCrystal_mfront_interface =
    "generic";

MGIS_SHAREDOBJ const char* FiniteStrainSingleCrystal_src =
    "FiniteStrainSingleCrystal.mfront";

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_nModellingHypotheses =
    1u;

MGIS_SHAREDOBJ const char* FiniteStrainSingleCrystal_ModellingHypotheses[1u] =
    {"Tridimensional"};

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_BehaviourType = 2u;

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_BehaviourKinematic =
    3u;

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_SymmetryType = 1u;

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_ElasticSymmetryType =
    1u;

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_OrthotropyManagementPolicy = 2u;

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_UsableInPurelyImplicitResolution =
        1;

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_nDrivingVariables = 1;
MGIS_SHAREDOBJ int FiniteStrainSingleCrystal_DrivingVariablesTypes[1] = {3};
MGIS_SHAREDOBJ const char * FiniteStrainSingleCrystal_DrivingVariables[1] = {"DeformationGradient"};

MGIS_SHAREDOBJ unsigned short FiniteStrainSingleCrystal_nThermodynamicForces = 1;
MGIS_SHAREDOBJ int FiniteStrainSingleCrystal_ThermodynamicForcesTypes[1] = {1};
MGIS_SHAREDOBJ const char * FiniteStrainSingleCrystal_ThermodynamicForces[1] = {"Stress"};

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nMaterialProperties = 7u;

MGIS_SHAREDOBJ const char*
    FiniteStrainSingleCrystal_Tridimensional_MaterialProperties[7u] = {
        "m", "K", "C", "R0", "Q", "b", "d1"};

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nInternalStateVariables = 37;
MGIS_SHAREDOBJ const char*
    FiniteStrainSingleCrystal_Tridimensional_InternalStateVariables[37] = {
        "g[0]",  "g[1]", "g[2]",  "g[3]",  "g[4]", "g[5]", "g[6]", "g[7]",
        "g[8]",  "g[9]", "g[10]", "g[11]", "Fe",   "p[0]", "p[1]", "p[2]",
        "p[3]",  "p[4]", "p[5]",  "p[6]",  "p[7]", "p[8]", "p[9]", "p[10]",
        "p[11]", "a[0]", "a[1]",  "a[2]",  "a[3]", "a[4]", "a[5]", "a[6]",
        "a[7]",  "a[8]", "a[9]",  "a[10]", "a[11]"};
MGIS_SHAREDOBJ int
    FiniteStrainSingleCrystal_Tridimensional_InternalStateVariablesTypes[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nExternalStateVariables = 0;
MGIS_SHAREDOBJ const char* const*
    FiniteStrainSingleCrystal_Tridimensional_ExternalStateVariables = nullptr;

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_nParameters = 12;
MGIS_SHAREDOBJ const char*
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
MGIS_SHAREDOBJ int
    FiniteStrainSingleCrystal_Tridimensional_ParametersTypes[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2};

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_theta_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_epsilon_ParameterDefaultValue =
        1e-11;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h1_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h2_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h3_ParameterDefaultValue = 0.6;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h4_ParameterDefaultValue = 12.3;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h5_ParameterDefaultValue = 1.6;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_h6_ParameterDefaultValue = 1.8;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_minimal_time_step_scaling_factor_ParameterDefaultValue =
        0.1;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_maximal_time_step_scaling_factor_ParameterDefaultValue =
        1.79769e+308;

MGIS_SHAREDOBJ double
    FiniteStrainSingleCrystal_Tridimensional_numerical_jacobian_epsilon_ParameterDefaultValue =
        1e-12;

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_iterMax_ParameterDefaultValue =
        100;

MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_requiresStiffnessTensor = 1;
MGIS_SHAREDOBJ unsigned short
    FiniteStrainSingleCrystal_Tridimensional_requiresThermalExpansionCoefficientTensor =
        0;

}  // end of extern "C"
