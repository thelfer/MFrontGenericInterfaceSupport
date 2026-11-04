macro(mgis_fortran_module file)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${file}"
    DESTINATION "${mgis_include_install_path}")
endmacro(mgis_fortran_module)

mgis_fortran_module(mgis_fortran_utilities.mod)
mgis_fortran_module(mgis.mod)
mgis_fortran_module(mgis_behaviour.mod)
mgis_fortran_module(mgis_model.mod)
