# preparing the environment
add_custom_target(mgis_python_bindings_checkfiles)
if(NOT CMAKE_CONFIGURATION_TYPES)
  add_custom_command(
  TARGET mgis_python_bindings_checkfiles POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E remove_directory mgis
    COMMAND ${CMAKE_COMMAND} -E make_directory mgis
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/bindings/python/pybind11/mgis/__init__.py" "mgis/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis>" "mgis/$<TARGET_FILE_NAME:py_mgis>"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis_behaviour>" "mgis/$<TARGET_FILE_NAME:py_mgis_behaviour>"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis_model>" "mgis/$<TARGET_FILE_NAME:py_mgis_model>"
    )
endif(NOT CMAKE_CONFIGURATION_TYPES)
add_dependencies(mgis_python_bindings_checkfiles py_mgis_behaviour)
add_dependencies(check mgis_python_bindings_checkfiles)

# defining the tests

function(test_python_bindings test)
  if(NOT CMAKE_CONFIGURATION_TYPES)
    add_test(NAME ${test}_py
      COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/${test}.py")
    set_property(TEST ${test}_py
      PROPERTY ENVIRONMENT
      PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python/pybind11/tests:$ENV{PYTHONPATH}
      MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>
      MGIS_TEST_MODELS_LIBRARY=$<TARGET_FILE:ModelTest>
      MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
  endif(NOT CMAKE_CONFIGURATION_TYPES)
endfunction(test_python_bindings)

test_python_bindings(MFrontGenericBehaviourInterfaceTest)
test_python_bindings(MFrontGenericBehaviourInterfaceTest2)
test_python_bindings(MFrontGenericBehaviourInterfaceTest3)
test_python_bindings(BoundsCheckTest)
test_python_bindings(ParameterTest)
test_python_bindings(IntegrateTest)
test_python_bindings(IntegrateTest2)
test_python_bindings(IntegrateTest3)
test_python_bindings(IntegrateTest4)
test_python_bindings(IntegrateTest5)
test_python_bindings(IntegrateTest6)
test_python_bindings(ExternalStateVariableTest)
test_python_bindings(InitializeFunctionTest)
test_python_bindings(PostProcessingTest)

