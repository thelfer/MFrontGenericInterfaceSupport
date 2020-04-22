/*!
* \file   StationaryHeatTransfer-generic.cxx
* \brief  This file implements the umat interface for the StationaryHeatTransfer behaviour law
* \author Thomas Helfer
* \date   15 / 02 / 2019
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#include<iostream>
#include<cstdlib>
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Material/StationaryHeatTransfer.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/StationaryHeatTransfer-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
StationaryHeatTransfer_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
StationaryHeatTransfer_build_id = "";

MFRONT_SHAREDOBJ const char* 
StationaryHeatTransfer_mfront_ept = "StationaryHeatTransfer";

MFRONT_SHAREDOBJ const char* 
StationaryHeatTransfer_tfel_version = "3.4.0-dev";

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
StationaryHeatTransfer_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
StationaryHeatTransfer_src = "StationaryHeatTransfer.mfront";

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
StationaryHeatTransfer_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nGradients = 1;

MFRONT_SHAREDOBJ int StationaryHeatTransfer_GradientsTypes[1] = {2};
MFRONT_SHAREDOBJ const char * StationaryHeatTransfer_Gradients[1] = {"TemperatureGradient"};
MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int StationaryHeatTransfer_ThermodynamicForcesTypes[1] = {2};
MFRONT_SHAREDOBJ const char * StationaryHeatTransfer_ThermodynamicForces[1] = {"HeatFlux"};
MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nTangentOperatorBlocks = 4;

MFRONT_SHAREDOBJ const char * StationaryHeatTransfer_TangentOperatorBlocks[4] = {"HeatFlux",
"TemperatureGradient","HeatFlux","Temperature"};
MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_BehaviourType = 0u;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_BehaviourKinematic = 0u;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_UsableInPurelyImplicitResolution = 0;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *StationaryHeatTransfer_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nInternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * StationaryHeatTransfer_InternalStateVariables = nullptr;

MFRONT_SHAREDOBJ const int * StationaryHeatTransfer_InternalStateVariablesTypes = nullptr;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * StationaryHeatTransfer_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_nParameters = 4;
MFRONT_SHAREDOBJ const char * StationaryHeatTransfer_Parameters[4] = {"A",
"B","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor"};
MFRONT_SHAREDOBJ int StationaryHeatTransfer_ParametersTypes [] = {0,0,0,0};

MFRONT_SHAREDOBJ double StationaryHeatTransfer_A_ParameterDefaultValue = 0.0375;

MFRONT_SHAREDOBJ double StationaryHeatTransfer_B_ParameterDefaultValue = 0.0002165;

MFRONT_SHAREDOBJ double StationaryHeatTransfer_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double StationaryHeatTransfer_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short StationaryHeatTransfer_requiresThermalExpansionCoefficientTensor = 0;

MFRONT_SHAREDOBJ void
StationaryHeatTransfer_setOutOfBoundsPolicy(const int p){
if(p==0){
StationaryHeatTransfer_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
StationaryHeatTransfer_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
StationaryHeatTransfer_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "StationaryHeatTransfer_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
StationaryHeatTransfer_setParameter(const char *const key,const double value){
using tfel::material::StationaryHeatTransferParametersInitializer;
auto& i = StationaryHeatTransferParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int StationaryHeatTransfer_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = StationaryHeatTransfer<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, StationaryHeatTransfer_getOutOfBoundsPolicy());
return r;
} // end of StationaryHeatTransfer_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int StationaryHeatTransfer_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = StationaryHeatTransfer<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, StationaryHeatTransfer_getOutOfBoundsPolicy());
return r;
} // end of StationaryHeatTransfer_Axisymmetrical

MFRONT_SHAREDOBJ int StationaryHeatTransfer_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = StationaryHeatTransfer<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, StationaryHeatTransfer_getOutOfBoundsPolicy());
return r;
} // end of StationaryHeatTransfer_PlaneStrain

MFRONT_SHAREDOBJ int StationaryHeatTransfer_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = StationaryHeatTransfer<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, StationaryHeatTransfer_getOutOfBoundsPolicy());
return r;
} // end of StationaryHeatTransfer_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int StationaryHeatTransfer_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr const auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = StationaryHeatTransfer<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, StationaryHeatTransfer_getOutOfBoundsPolicy());
return r;
} // end of StationaryHeatTransfer_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

