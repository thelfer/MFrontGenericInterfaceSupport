add_library(py_mgis_behaviour MODULE
 behaviour-module.cxx
 Hypothesis.cxx
 Variable.cxx
 Behaviour.cxx
 State.cxx
 BehaviourData.cxx
 Integrate.cxx)

if(WIN32)
  set_target_properties(py_mgis_behaviour PROPERTIES
    COMPILE_FLAGS "-DHAVE_ROUND")
  set_target_properties(py_mgis_behaviour PROPERTIES SUFFIX ".pyd")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(py_mgis_behaviour PROPERTIES SUFFIX ".so")
endif(WIN32)
set_target_properties(py_mgis_behaviour PROPERTIES PREFIX "")
set_target_properties(py_mgis_behaviour PROPERTIES OUTPUT_NAME behaviour)

# include directories
target_include_directories(py_mgis_behaviour
  SYSTEM
  PRIVATE "${CMAKE_SOURCE_DIR}/bindings/python/include"
  PRIVATE "${Boost_INCLUDE_DIRS}"
  PRIVATE "${PYTHON_INCLUDE_DIRS}")

# link flags
if(MGIS_USES_CONAN)
  target_link_libraries(py_mgis_behaviour
    PRIVATE
    MFrontGenericInterface
    ${CONAN_LIBS} ${PYTHON_LIBRARIES})
else(MGIS_USES_CONAN)
  target_link_libraries(py_mgis_behaviour
    PRIVATE
    MFrontGenericInterface
    ${Boost_NUMPY_LIBRARY} ${Boost_PYTHON_LIBRARY} ${PYTHON_LIBRARIES})
endif(MGIS_USES_CONAN)

if(MGIS_APPEND_SUFFIX)
  if(WIN32)
    install(TARGETS py_mgis_behaviour
      DESTINATION bin/${PYTHON_LIBRARY}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES}
      COMPONENT python_bindings)
  else(WIN32)
    install(TARGETS py_mgis_behaviour
      DESTINATION lib${LIB_SUFFIX}/${PYTHON_LIBRARY}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES}
      COMPONENT python_bindings)
  endif(WIN32)
else(MGIS_APPEND_SUFFIX)
  if(WIN32)
    install(TARGETS py_mgis_behaviour
      DESTINATION bin/${PYTHON_LIBRARY}/site-packages/mgis
      COMPONENT python_bindings)
  else(WIN32)
    install(TARGETS py_mgis_behaviour
      DESTINATION lib${LIB_SUFFIX}/${PYTHON_LIBRARY}/site-packages/mgis
      COMPONENT python_bindings)
  endif(WIN32)
endif(MGIS_APPEND_SUFFIX)
