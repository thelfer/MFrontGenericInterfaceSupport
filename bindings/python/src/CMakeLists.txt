function(mgis_python_module name output_name)
  if(${ARGC} LESS 3)
    message(FATAL_ERROR "mgis_python_module: no source specified")
  endif(${ARGC} LESS 3)
  set(module "py_${name}")
  add_library(${module} MODULE ${ARGN})
  if(enable-broken-boost-python-module-visibility-handling)
    target_compile_definitions(${module}
      PRIVATE BOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY)
  endif(enable-broken-boost-python-module-visibility-handling)
  if(WIN32)
    set_target_properties(${module} PROPERTIES
      COMPILE_FLAGS "-DHAVE_ROUND")
    set_target_properties(${module} PROPERTIES SUFFIX ".pyd")
  elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${module} PROPERTIES SUFFIX ".so")
  endif(WIN32)
  set_target_properties(${module} PROPERTIES PREFIX "")
  set_target_properties(${module} PROPERTIES OUTPUT_NAME "${output_name}")
  # include directories
  target_include_directories(${module}
    SYSTEM
    PRIVATE "${CMAKE_SOURCE_DIR}/bindings/python/include"
    PRIVATE "${Boost_INCLUDE_DIRS}"
    PRIVATE "${PYTHON_INCLUDE_DIRS}"
    PRIVATE "${PYTHON_NUMPY_INCLUDE_DIR}")
  # link flags
  if(MGIS_USES_CONAN)
    target_link_libraries(${module}
      PRIVATE
      MFrontGenericInterface
      ${CONAN_LIBS} ${PYTHON_LIBRARIES})
  else(MGIS_USES_CONAN)
    target_link_libraries(${module}
      PRIVATE
      MFrontGenericInterface
      ${Boost_NUMPY_LIBRARY} ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARIES})
  endif(MGIS_USES_CONAN)
  # installation
  install(TARGETS ${module}
    DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
    COMPONENT python_bindings)
endfunction(mgis_python_module)

mgis_python_module(mgis _mgis
 mgis-module.cxx
 ThreadPool.cxx)

mgis_python_module(mgis_behaviour behaviour
 behaviour-module.cxx
 NumPySupport.cxx
 Hypothesis.cxx
 Variable.cxx
 Behaviour.cxx
 State.cxx
 BehaviourData.cxx
 BehaviourDataView.cxx
 MaterialDataManager.cxx
 MaterialStateManager.cxx
 Integrate.cxx
 FiniteStrainSupport.cxx)
