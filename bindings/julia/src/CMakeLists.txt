# store the location of the julia modules in the MGIS_JULIA_MODULES_INSTALL_DIRECTORY variable
if(NOT DEFINED MGIS_JULIA_MODULES_INSTALL_DIRECTORY)
  if(MGIS_APPEND_SUFFIX)
    if(WIN32)
      set(MGIS_JULIA_MODULES_INSTALL_DIRECTORY
          bin/julia/mgis_${MGIS_SUFFIX_FOR_JULIA_MODULES})
    else(WIN32)
      set(MGIS_JULIA_MODULES_INSTALL_DIRECTORY
          "lib${LIB_SUFFIX}/julia/mgis_${MGIS_SUFFIX_FOR_JULIA_MODULES}")
    endif(WIN32)
  else(MGIS_APPEND_SUFFIX)
    if(WIN32)
      set(MGIS_JULIA_MODULES_INSTALL_DIRECTORY
          "bin/julia/mgis")
    else(WIN32)
      set(MGIS_JULIA_MODULES_INSTALL_DIRECTORY
          "lib${LIB_SUFFIX}/julia/mgis")
    endif(WIN32)
  endif(MGIS_APPEND_SUFFIX)
endif(NOT DEFINED MGIS_JULIA_MODULES_INSTALL_DIRECTORY)

function(mgis_julia_module name)
  if(${ARGC} LESS 2)
    message(FATAL_ERROR "mgis_julia_module: no source specified")
  endif(${ARGC} LESS 2)
  set(module "${name}-julia")
  add_library(${module} MODULE ${ARGN})
  set_target_properties(${module} PROPERTIES PREFIX "")
  # include directories
  target_include_directories(${module}
    PRIVATE "${PROJECT_SOURCE_DIR}/bindings/julia/include")
  # link flags
  if(MGIS_USES_CONAN)
    target_link_libraries(${module}
      PUBLIC MFrontGenericInterface
      PRIVATE JlCxx::cxxwrap_julia
      ${CONAN_LIBS})
  else(MGIS_USES_CONAN)
    target_link_libraries(${module}
      PUBLIC MFrontGenericInterface
      PRIVATE JlCxx::cxxwrap_julia)
  endif(MGIS_USES_CONAN)
  # installation
  install(TARGETS ${module}
    DESTINATION ${MGIS_JULIA_MODULES_INSTALL_DIRECTORY}
    COMPONENT julia_bindings)
endfunction(mgis_julia_module)

mgis_julia_module(mgis
 mgis-module.cxx
 ThreadPool.cxx
 behaviour-module.cxx
 Hypothesis.cxx
 Variable.cxx
 Behaviour.cxx
 State.cxx
 BehaviourData.cxx
 BehaviourDataView.cxx
 MaterialStateManager.cxx
 MaterialDataManager.cxx
 Integrate.cxx
)

