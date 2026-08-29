string(REPLACE "\\;" ";" MGIS_PATH_STRING "$ENV{PATH}")
string(REPLACE ";" "\\;" MGIS_PATH_STRING "${MGIS_PATH_STRING}")

# =============================================================================
# Helper Macros
# =============================================================================

# Set Windows-specific test environment
macro(set_test_environment test_name)
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test_name} PROPERTY
      ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
  endif()
endmacro()

# Add compiler flags for function tests (GNU < 14.2.0 and NVHPC)
macro(add_function_test_flags name)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
      target_compile_definitions(${name} PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(${name} PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif()
endmacro()

# Add a simple test executable
macro(add_simple_test name source)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface)
  add_test(NAME ${name} COMMAND ${name})
  add_dependencies(build-tests ${name})
  set_test_environment(${name})
endmacro()

# Add a test that depends on a library
macro(add_library_test name source lib)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface)
  add_test(NAME ${name} COMMAND ${name} "$<TARGET_FILE:${lib}>")
  add_dependencies(build-tests ${name})
  set_property(TEST ${name} PROPERTY DEPENDS ${lib})
  set_test_environment(${name})
endmacro()

# Add a library test with TFEL_VERSION definition
macro(add_library_test_with_tfel_version name source lib)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface)
  target_compile_definitions(${name} PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
  add_test(NAME ${name} COMMAND ${name} "$<TARGET_FILE:${lib}>")
  add_dependencies(build-tests ${name})
  set_property(TEST ${name} PROPERTY DEPENDS ${lib})
  set_test_environment(${name})
endmacro()

# Add a test with specific behaviour
macro(add_behaviour_test name source lib behaviour)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface)
  add_test(NAME ${name} COMMAND ${name} "$<TARGET_FILE:${lib}>" "${behaviour}")
  add_dependencies(build-tests ${name})
  set_property(TEST ${name} PROPERTY DEPENDS ${lib})
  set_test_environment(${name})
endmacro()

# Add a test with TFEL_VERSION definition
macro(add_tfel_version_test name source lib behaviour)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface)
  target_compile_definitions(${name} PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
  add_test(NAME ${name} COMMAND ${name} "$<TARGET_FILE:${lib}>" "${behaviour}")
  add_dependencies(build-tests ${name})
  set_property(TEST ${name} PROPERTY DEPENDS ${lib})
  set_test_environment(${name})
endmacro()

# Add a test with TFEL tests library
macro(add_tfel_test name source)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface tfel::TFELTests)
  add_test(NAME ${name} COMMAND ${name})
  add_dependencies(build-tests ${name})
  set_test_environment(${name})
endmacro()

# Add a test with TFEL tests library and behaviour library
macro(add_tfel_library_test name source lib)
  add_executable(${name} EXCLUDE_FROM_ALL ${source})
  target_link_libraries(${name} PRIVATE MFrontGenericInterface tfel::TFELTests)
  add_test(NAME ${name} COMMAND ${name} "$<TARGET_FILE:${lib}>")
  add_dependencies(build-tests ${name})
  set_property(TEST ${name} PROPERTY DEPENDS ${lib})
  set_test_environment(${name})
endmacro()

# ===========================================================================
# Simple Tests
# ===========================================================================

add_simple_test(IntegrationFailureDebugTest IntegrationFailureDebugTest.cxx)
add_simple_test(InvalidResultTest InvalidResultTest.cxx)
add_simple_test(ErrorBacktraceTest ErrorBacktraceTest.cxx)

# =============================================================================
# MFront Behaviour Libraries (unchanged from original)
# =============================================================================

if(MGIS_HAVE_MFRONT_SUPPORT)

  mfront_material_properties_check_library(
    MaterialPropertyTest
    MaterialPropertyBoundCheck
    MaterialPropertyBoundCheck2
    MaterialPropertyBoundCheck3
    MaterialPropertyErrnoCheck
    MaterialPropertyErrnoCheck2
    Inconel600_YoungModulus)

  mfront_behaviours_check_library(
    BehaviourTest
    Elasticity
    OrthotropicElasticity
    Norton
    NortonFailure
    Plasticity
    Gurson
    BoundsCheckTest
    ParameterTest
    FiniteStrainSingleCrystal
    StandardElastoViscoPlasticityPlasticityTest11
    TensorialExternalStateVariableTest
    InitializeFunctionTest
    PostProcessingTest
    ComputeSpeedOfSound)

  mfront_behaviours_check_library(ModelTest ode_rk54)

  # Tests that depend on behaviour libraries
  add_library_test_with_tfel_version(LoadMaterialPropertyTest LoadMaterialPropertyTest.cxx MaterialPropertyTest)

  # Tests with TFEL_VERSION
  add_tfel_version_test(MFrontGenericBehaviourInterfaceTest
    MFrontGenericBehaviourInterfaceTest.cxx BehaviourTest Gurson)
  add_tfel_version_test(MFrontGenericBehaviourInterfaceTest2
    MFrontGenericBehaviourInterfaceTest2.cxx BehaviourTest FiniteStrainSingleCrystal)
  add_tfel_version_test(MFrontGenericBehaviourInterfaceTest3
    MFrontGenericBehaviourInterfaceTest3.cxx BehaviourTest FiniteStrainSingleCrystal)
  add_tfel_version_test(BoundsCheckTest BoundsCheckTest.cxx BehaviourTest BoundsCheckTest)
  add_tfel_version_test(ParameterTest ParameterTest.cxx BehaviourTest ParameterTest)

  # Other behaviour tests
  add_library_test(UpdatePolicyTest UpdatePolicyTest.cxx BehaviourTest)
  add_library_test(IntegrateTest IntegrateTest.cxx BehaviourTest)
  add_library_test(IntegrateTest2 IntegrateTest2.cxx BehaviourTest)
  add_library_test(IntegrateTest2b IntegrateTest2b.cxx BehaviourTest)
  add_library_test(IntegrateTest3 IntegrateTest3.cxx BehaviourTest)
  add_library_test(IntegrateTest3b IntegrateTest3b.cxx BehaviourTest)
  add_library_test(IntegrateTest4 IntegrateTest4.cxx ModelTest)
  add_library_test(IntegrateTest5 IntegrateTest5.cxx ModelTest)
  add_library_test(RotateFunctionsTest RotateFunctionsTest.cxx BehaviourTest)
  add_library_test(IntegrateDebugTest IntegrateDebugTest.cxx BehaviourTest)
  add_tfel_version_test(ExternalStateVariableTest ExternalStateVariableTest.cxx
    BehaviourTest TensorialExternalStateVariableTest)
  add_tfel_version_test(InitializeFunctionTest InitializeFunctionTest.cxx
    BehaviourTest InitializeFunctionTest)
  add_tfel_version_test(PostProcessingTest PostProcessingTest.cxx
    BehaviourTest PostProcessingTest)
  add_library_test(ComputeSpeedOfSoundTest ComputeSpeedOfSoundTest.cxx BehaviourTest)
  add_library_test(ComputeSpeedOfSoundTest2 ComputeSpeedOfSoundTest2.cxx BehaviourTest)

endif(MGIS_HAVE_MFRONT_SUPPORT)

# ===========================================================================
# TFEL-dependent Tests
# ===========================================================================

if(MGIS_HAVE_TFEL)
  add_tfel_test(ContextTest ContextTest.cxx)
  add_tfel_test(ProfilingTest ProfilingTest.cxx)
  add_tfel_library_test(OptionalReferenceTest OptionalReferenceTest.cxx BehaviourTest)
  add_tfel_library_test(ConstructTest ConstructTest.cxx BehaviourTest)
  add_tfel_library_test(InvokeTest InvokeTest.cxx BehaviourTest)

  if(enable-mgis-function)
    add_tfel_test(SharedSpaceTest SharedSpaceTest.cxx)

    # Function tests with compiler-specific flags
    add_tfel_test(FunctionTest FunctionTest.cxx)
    add_function_test_flags(FunctionTest)

    add_tfel_test(EvaluatorsTest EvaluatorsTest.cxx)
    add_function_test_flags(EvaluatorsTest)

    add_tfel_test(UniformEvaluatorsTest UniformEvaluatorsTest.cxx)
    add_function_test_flags(UniformEvaluatorsTest)

    add_tfel_test(CoalescedMemoryAccessFunctionViewTest CoalescedMemoryAccessFunctionViewTest.cxx)
    add_function_test_flags(CoalescedMemoryAccessFunctionViewTest)

    add_tfel_test(StridedCoalescedMemoryAccessFunctionViewTest StridedCoalescedMemoryAccessFunctionViewTest.cxx)
    add_function_test_flags(StridedCoalescedMemoryAccessFunctionViewTest)

    add_tfel_test(MechanicalEvaluatorsTest MechanicalEvaluatorsTest.cxx)
    add_function_test_flags(MechanicalEvaluatorsTest)

    add_tfel_test(QuantityFunctionTest QuantityFunctionTest.cxx)
    add_function_test_flags(QuantityFunctionTest)

    add_tfel_test(TensorialFunctionTest TensorialFunctionTest.cxx)
    add_function_test_flags(TensorialFunctionTest)

    add_tfel_test(ReduceTest ReduceTest.cxx)

    # MaterialFunctionManagerTest - executable created here, test registered conditionally
    add_executable(MaterialFunctionManagerTest EXCLUDE_FROM_ALL MaterialFunctionManagerTest.cxx)
    target_link_libraries(MaterialFunctionManagerTest PRIVATE MFrontGenericInterface tfel::TFELTests)
    add_dependencies(build-tests MaterialFunctionManagerTest)
    set_test_environment(MaterialFunctionManagerTest)

    if(enable-cuda-algorithms)
      add_executable(FunctionTest-cuda EXCLUDE_FROM_ALL FunctionTest-cuda.cu)
      target_compile_options(FunctionTest-cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${MGIS_CUDA_FLAGS}>)
      target_link_libraries(FunctionTest-cuda PRIVATE MFrontGenericInterface tfel::TFELTests)
      add_test(NAME FunctionTest-cuda COMMAND FunctionTest-cuda "$<TARGET_FILE:BehaviourTest>")
      add_dependencies(build-tests FunctionTest-cuda)
      set_property(TEST FunctionTest-cuda PROPERTY DEPENDS BehaviourTest)
      set_test_environment(FunctionTest-cuda)
    endif()

    if(MGIS_HAVE_MFRONT_SUPPORT)
      add_test(NAME MaterialFunctionManagerTest
        COMMAND MaterialFunctionManagerTest "$<TARGET_FILE:BehaviourTest>")
      set_property(TEST MaterialFunctionManagerTest PROPERTY DEPENDS BehaviourTest)
    endif()
  endif(enable-mgis-function)

  if(MGIS_HDF5_SUPPORT)
    add_library_test(HDF5Test HDF5Test.cxx BehaviourTest)
  endif()

endif(MGIS_HAVE_TFEL)
