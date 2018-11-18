macro(mgis_fortran_module file)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${file}"
    DESTINATION "include")
endmacro(mgis_fortran_module)

mgis_fortran_module(mgis.mod)
mgis_fortran_module(mgis_behaviour.mod)
