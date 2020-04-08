function(mgis_fenics_python_script file)
  install(FILES ${file}
          DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}/fencis
          COMPONENT python_bindings)
endfunction(mgis_fenics_python_script)

mgis_fenics_python_script(__init__.py)
