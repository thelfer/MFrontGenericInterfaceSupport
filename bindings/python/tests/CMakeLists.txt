# preparing the environment
add_custom_target(mgis_python_bindings_checkfiles)
if(CMAKE_CONFIGURATION_TYPES)
  # add_custom_command(
  # TARGET mgis_python_bindings_checkfiles
  #   foreach(conf ${CMAKE_CONFIGURATION_TYPES})
  #     COMMAND ${CMAKE_COMMAND} -E remove_directory "${conf}/mgis"
  #     COMMAND ${CMAKE_COMMAND} -E make_directory "${conf}/mgis"
  #     COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/bindings/python/src/__init__.py" "${conf}/mgis/__init__.py"
  #     COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis_behaviour>" "${conf}/mgis/$<TARGET_FILE_NAME:py_mgis_behaviour>"
  #   endforeach(conf ${CMAKE_CONFIGURATION_TYPES})
  # )
else(CMAKE_CONFIGURATION_TYPES)
  add_custom_command(
  TARGET mgis_python_bindings_checkfiles
    COMMAND ${CMAKE_COMMAND} -E remove_directory mgis
    COMMAND ${CMAKE_COMMAND} -E make_directory mgis
    COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/bindings/python/src/__init__.py" "mgis/__init__.py"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis>" "mgis/$<TARGET_FILE_NAME:py_mgis>"
    COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis_behaviour>" "mgis/$<TARGET_FILE_NAME:py_mgis_behaviour>"
    )
endif(CMAKE_CONFIGURATION_TYPES)
add_dependencies(mgis_python_bindings_checkfiles py_mgis_behaviour)
add_dependencies(check mgis_python_bindings_checkfiles)

# defining the tests

function(test_python_bindings test)
  if(CMAKE_CONFIGURATION_TYPES)
    foreach(conf ${CMAKE_CONFIGURATION_TYPES})
      add_test(NAME ${test}_${conf}_py
	COMMAND ${PYTHON_EXECUTABLE} ${test}.py "$<TARGET_FILE:BehaviourTest>"
	CONFIGURATION ${conf})
      set_tests_properties(${test}_${conf}_py
    	PROPERTIES ENVIRONMENT
        PYTHONPATH="${PROJECT_BINARY_DIR}/bindings/python/tests/;$ENV{PYTHONPATH}"
        MGIS_TEST_BEHAVIOURS_LIBRARY="$<TARGET_FILE:BehaviourTest>")
    endforeach(conf ${CMAKE_CONFIGURATION_TYPES})
  else(CMAKE_CONFIGURATION_TYPES)
    add_test(NAME ${test}_py
      COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/${test}.py")
    set_property(TEST ${test}_py
      PROPERTY ENVIRONMENT
      PYTHONPATH=${PROJECT_BINARY_DIR}/bindings/python/tests/;$ENV{PYTHONPATH}
      MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>
      MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
  endif(CMAKE_CONFIGURATION_TYPES)
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

