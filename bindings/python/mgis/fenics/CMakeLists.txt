function(mgis_fenics_python_script file)
  install(FILES ${file}
          DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}/fenics
          COMPONENT python_bindings)
endfunction(mgis_fenics_python_script)

mgis_fenics_python_script(__init__.py)
mgis_fenics_python_script(gradient_flux.py)
mgis_fenics_python_script(utils.py)
mgis_fenics_python_script(nonlinear_material.py)
mgis_fenics_python_script(nonlinear_problem.py)
