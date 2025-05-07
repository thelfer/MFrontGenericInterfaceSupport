mgis_library(MFrontGenericInterface
  Raise.cxx
  VerbosityLevel.cxx
  ErrorBacktrace.cxx
  Context.cxx
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
  BehaviourData.cxx
  MaterialStateManager.cxx
  MaterialDataManager.cxx
  Integrate.cxx
  FiniteStrainSupport.cxx
  Model.cxx
  BasicLinearSpace.cxx
  BasicLinearQuadratureSpace.cxx
  Function.cxx
  Evaluators.cxx
  MechanicalEvaluators.cxx
  )

target_include_directories(MFrontGenericInterface
   PUBLIC 
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:include>)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_DEBUG)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if(enable-static)
  target_include_directories(MFrontGenericInterface-static
    PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(MFrontGenericInterface-static PUBLIC MGIS_DEBUG)
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
endif(enable-static)

if(MGIS_HAVE_TFEL)
  target_compile_definitions(MFrontGenericInterface PUBLIC MGIS_HAVE_TFEL)
  target_link_libraries(MFrontGenericInterface PUBLIC tfel::TFELMath)
endif(MGIS_HAVE_TFEL)

if(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY} Threads::Threads)
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE ${MGIS_DL_LIBRARY} Threads::Threads)
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
      PRIVATE Threads::Threads)
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE Threads::Threads)
    endif(enable-static)
  endif(Threads_FOUND)
endif(UNIX)
