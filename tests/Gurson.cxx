/*!
 * \file   Gurson.cxx
 * \brief
 * \author Thomas Helfer
 * \date   20/06/2018
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

MGIS_SHAREDOBJ const char *gurson_mfront_ept = "gurson";

MGIS_SHAREDOBJ const char *gurson_tfel_version = "3.2.0-dev";

MGIS_SHAREDOBJ unsigned short gurson_mfront_mkt = 1u;

MGIS_SHAREDOBJ const char *gurson_src = "Gurson.mfront";

MGIS_SHAREDOBJ const char *gurson_mfront_interface = "generic";

MGIS_SHAREDOBJ unsigned short gurson_nModellingHypotheses = 4u;

MGIS_SHAREDOBJ const char *gurson_ModellingHypotheses[4u] = {
    "Axisymmetrical", "PlaneStrain", "GeneralisedPlaneStrain",
    "Tridimensional"};

MGIS_SHAREDOBJ unsigned short gurson_BehaviourType = 1u;

MGIS_SHAREDOBJ unsigned short gurson_BehaviourKinematic = 1u;

MGIS_SHAREDOBJ unsigned short gurson_SymmetryType = 0u;

MGIS_SHAREDOBJ unsigned short gurson_ElasticSymmetryType = 0u;

MGIS_SHAREDOBJ unsigned short gurson_savesTangentOperator = 0;
MGIS_SHAREDOBJ unsigned short gurson_UsableInPurelyImplicitResolution = 1;

MGIS_SHAREDOBJ unsigned short gurson_nMaterialProperties = 0u;

MGIS_SHAREDOBJ const char *const *gurson_MaterialProperties = nullptr;

MGIS_SHAREDOBJ unsigned short gurson_nInternalStateVariables = 4;
MGIS_SHAREDOBJ const char *gurson_InternalStateVariables[4] = {
    "ElasticStrain", "EquivalentPlasticStrain", "MatrixEquivalentPlasticStrain",
    "Porosity"};
MGIS_SHAREDOBJ int gurson_InternalStateVariablesTypes[] = {1, 0, 0, 0};

MGIS_SHAREDOBJ unsigned short gurson_nExternalStateVariables = 0;
MGIS_SHAREDOBJ const char *const *gurson_ExternalStateVariables = nullptr;

MGIS_SHAREDOBJ unsigned short gurson_nParameters = 16;
MGIS_SHAREDOBJ const char *gurson_Parameters[16] = {
    "epsilon",
    "theta",
    "FirstGursonParameter",
    "SecondGursonParameter",
    "ThirdGursonParameter",
    "ffb",
    "fmax",
    "fc",
    "sm",
    "fNH",
    "en",
    "sn",
    "minimal_time_step_scaling_factor",
    "maximal_time_step_scaling_factor",
    "numerical_jacobian_epsilon",
    "iterMax"};
MGIS_SHAREDOBJ int gurson_ParametersTypes[] = {0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 2};

MGIS_SHAREDOBJ double gurson_epsilon_ParameterDefaultValue = 1e-12;

MGIS_SHAREDOBJ double gurson_theta_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double gurson_FirstGursonParameter_ParameterDefaultValue = 1.5;

MGIS_SHAREDOBJ double gurson_SecondGursonParameter_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double gurson_ThirdGursonParameter_ParameterDefaultValue = 1;

MGIS_SHAREDOBJ double gurson_ffb_ParameterDefaultValue = 0.5;

MGIS_SHAREDOBJ double gurson_fmax_ParameterDefaultValue = 0.18;

MGIS_SHAREDOBJ double gurson_fc_ParameterDefaultValue = 0.1;

MGIS_SHAREDOBJ double gurson_sm_ParameterDefaultValue = 4e+08;

MGIS_SHAREDOBJ double gurson_fNH_ParameterDefaultValue = 0.1;

MGIS_SHAREDOBJ double gurson_en_ParameterDefaultValue = 0.3;

MGIS_SHAREDOBJ double gurson_sn_ParameterDefaultValue = 0.1;

MGIS_SHAREDOBJ double
    gurson_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MGIS_SHAREDOBJ double
    gurson_maximal_time_step_scaling_factor_ParameterDefaultValue =
        1.79769e+308;

MGIS_SHAREDOBJ double gurson_numerical_jacobian_epsilon_ParameterDefaultValue =
    1e-13;

MGIS_SHAREDOBJ unsigned short gurson_iterMax_ParameterDefaultValue = 100;

MGIS_SHAREDOBJ long double gurson_EquivalentPlasticStrain_LowerPhysicalBound =
    0;

MGIS_SHAREDOBJ long double
    gurson_MatrixEquivalentPlasticStrain_LowerPhysicalBound = 0;

MGIS_SHAREDOBJ long double gurson_Porosity_LowerPhysicalBound = 0;

MGIS_SHAREDOBJ long double gurson_Porosity_UpperPhysicalBound = 1;

MGIS_SHAREDOBJ unsigned short gurson_requiresStiffnessTensor = 0;
MGIS_SHAREDOBJ unsigned short gurson_requiresThermalExpansionCoefficientTensor =
    0;

}  // end of extern "C"
