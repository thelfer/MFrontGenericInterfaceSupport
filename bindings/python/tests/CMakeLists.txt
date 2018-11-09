add_custom_target(mgis_python_bindings_checkfiles)
add_custom_command(
  TARGET mgis_python_bindings_checkfiles
  COMMAND ${CMAKE_COMMAND} -E remove_directory mgis
  COMMAND ${CMAKE_COMMAND} -E make_directory mgis
  COMMAND ${CMAKE_COMMAND} -E create_symlink "${PROJECT_SOURCE_DIR}/bindings/python/src/__init__.py" "mgis/__init__.py"
  COMMAND ${CMAKE_COMMAND} -E create_symlink "$<TARGET_FILE:py_mgis_behaviour>" "mgis/$<TARGET_FILE_NAME:py_mgis_behaviour>"
  )

add_dependencies(mgis_python_bindings_checkfiles py_mgis_behaviour)
add_dependencies(check mgis_python_bindings_checkfiles)
