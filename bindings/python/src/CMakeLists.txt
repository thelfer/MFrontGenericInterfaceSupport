# store the location of the python modules in the MGIS_PYTHON_MODULES_INSTALL_DIRECTORY variable
if(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        bin/${PYTHON_LIBRARY}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES})
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/${PYTHON_LIBRARY}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES}")
  endif(WIN32)
else(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "bin/${PYTHON_LIBRARY}/site-packages/mgis")
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/${PYTHON_LIBRARY}/site-packages/mgis")
  endif(WIN32)
endif(MGIS_APPEND_SUFFIX)

function(mgis_python_module name output_name)
  if(${ARGC} LESS 3)
    message(FATAL_ERROR "mgis_python_module: no source specified")
  endif(${ARGC} LESS 3)
  set(module "py_${name}")
  add_library(${module} MODULE ${ARGN})
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
    PRIVATE "${PYTHON_INCLUDE_DIRS}")
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
 Integrate.cxx)

install(FILES __init__.py
        DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
        COMPONENT python_bindings)

