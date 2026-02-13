set(MFrontGenericInterface_SOURCES)
list(APPEND MFrontGenericInterface_SOURCES
  Raise.cxx
  VerbosityLevel.cxx
  ErrorBacktrace.cxx
  Context.cxx
  Contract.cxx
  ThreadPool.cxx
  ThreadedTaskResult.cxx
  LibrariesManager.cxx
  Markdown.cxx
  MaterialProperty.cxx
  RotationMatrix.cxx
  Variable.cxx
  Hypothesis.cxx
  BehaviourDescription.cxx
  Behaviour.cxx
  State.cxx
  StateView.cxx
  BehaviourData.cxx
  BehaviourDataView.cxx
  MaterialStateManager.cxx
  MaterialDataManager.cxx
  Integrate.cxx
  FiniteStrainSupport.cxx
  Model.cxx
  BehaviourIntegrationFailureAnalyser.cxx)

if(MGIS_HAVE_TFEL)
  list(APPEND MFrontGenericInterface_SOURCES
       Database.cxx)
endif(MGIS_HAVE_TFEL)

if(enable-mgis-function)
  list(APPEND MFrontGenericInterface_SOURCES
    BasicLinearSpace.cxx
    BasicLinearQuadratureSpace.cxx
    FixedSizeView.cxx
    Function.cxx
    Evaluator.cxx
    UniformEvaluator.cxx
    CoalescedMemoryAccessFunctionView.cxx
    StridedCoalescedMemoryAccessFunctionView.cxx
    MaterialFunctionManager.cxx)
  if(MGIS_HAVE_TFEL)
    list(APPEND MFrontGenericInterface_SOURCES
         TensorView.cxx
         TensorModifier.cxx)
    if(enable-mgis-function-precompiled-operations)
      list(APPEND MFrontGenericInterface_SOURCES
           TensorOperations.cxx
           MechanicalOperations.cxx)
    endif(enable-mgis-function-precompiled-operations)
  endif(MGIS_HAVE_TFEL)
endif(enable-mgis-function)

if(MGIS_HDF5_SUPPORT)
  list(APPEND MFrontGenericInterface_SOURCES HDF5Support.cxx)
endif(MGIS_HDF5_SUPPORT)


mgis_library(MFrontGenericInterface
  ${MFrontGenericInterface_SOURCES})

target_include_directories(MFrontGenericInterface
   PUBLIC 
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:include>)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_DEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
if(enable-mgis-function)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_FUNCTION_SUPPORT)
endif(enable-mgis-function)
if(enable-exceptions)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING)
endif(enable-exceptions)
if(enable-parallel-stl-algorithms)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_USE_STL_PARALLEL_ALGORITHMS)
endif(enable-parallel-stl-algorithms)

if(enable-static)
  target_include_directories(MFrontGenericInterface-static
    PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_DEBUG)
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if(enable-mgis-function)
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_FUNCTION_SUPPORT)
  endif(enable-mgis-function)
  if(enable-parallel-stl-algorithms)
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_USE_STL_PARALLEL_ALGORITHMS)
  endif(enable-parallel-stl-algorithms)
  if(MGIS_HAVE_TFEL)
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_HAVE_TFEL)
  endif(MGIS_HAVE_TFEL)
endif(enable-static)

if(MGIS_HAVE_TFEL)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_HAVE_TFEL)
  target_link_libraries(MFrontGenericInterface PUBLIC tfel::TFELMath tfel::TFELMFrontDatabase)
  if(enable-static)                                                                                                                                    
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_HAVE_TFEL)
    if(TARGET tfel::TFELMath-static)                                                                
      target_link_libraries(MFrontGenericInterface-static PUBLIC tfel::TFELMath-static)
    else()
      target_link_libraries(MFrontGenericInterface-static PUBLIC tfel::TFELMath)
    endif()
    if(TARGET tfel::TFELMFrontDatabase-static)                                                                
      target_link_libraries(MFrontGenericInterface-static PUBLIC tfel::TFELMFrontDatabase-static)                                                
    else()
      target_link_libraries(MFrontGenericInterface-static PUBLIC tfel::TFELMFrontDatabase)                                              
    endif()
  endif(enable-static) 
endif(MGIS_HAVE_TFEL)

if(enable-openmp)
  target_link_libraries(MFrontGenericInterface PUBLIC OpenMP::OpenMP_CXX)
endif(enable-openmp)

if(MGIS_ADDITIONAL_LIBRARIES)
  target_link_libraries(MFrontGenericInterface PUBLIC ${MGIS_ADDITIONAL_LIBRARIES})
endif(MGIS_ADDITIONAL_LIBRARIES)

if(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PUBLIC  Threads::Threads
      PRIVATE ${MGIS_DL_LIBRARY})
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PUBLIC  Threads::Threads
        PRIVATE ${MGIS_DL_LIBRARY})
    endif(enable-static)
  else(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY})
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE ${MGIS_DL_LIBRARY})
    endif(enable-static)
  endif(Threads_FOUND)
else(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PUBLIC Threads::Threads)
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PUBLIC Threads::Threads)
    endif(enable-static)
  endif(Threads_FOUND)
endif(UNIX)

if(MGIS_HDF5_SUPPORT)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_HAVE_HDF5)
  target_compile_definitions(MFrontGenericInterface PUBLIC ${HDF5_DEFINITIONS})
  target_include_directories(MFrontGenericInterface SYSTEM  PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(MFrontGenericInterface PUBLIC ${HDF5_LIBRARIES})
  # target_link_libraries(MFrontGenericInterface PUBLIC hdf5::hdf5_cpp)
  if(enable-static)
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_HAVE_HDF5)
    target_compile_definitions(MFrontGenericInterface-static PUBLIC ${HDF5_DEFINITIONS})
    target_include_directories(MFrontGenericInterface-static SYSTEM  PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(MFrontGenericInterface-static PUBLIC ${HDF5_LIBRARIES})
    # target_link_libraries(MFrontGenericInterface-static PUBLIC hdf5::hdf5_cpp hdf5::hdf5)
  endif(enable-static)
endif(MGIS_HDF5_SUPPORT)

if(MGIS_ADDITIONAL_LINK_FLAGS)
  target_link_options(MFrontGenericInterface PUBLIC "${MGIS_ADDITIONAL_LINK_FLAGS}")
endif(MGIS_ADDITIONAL_LINK_FLAGS)
