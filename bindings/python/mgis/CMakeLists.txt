add_subdirectory(fenics)

if(MGIS_HAVE_TFEL AND TFEL_PYTHON_BINDINGS)
  install(FILES __init__-have-tfel.py
          DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
          COMPONENT python_bindings
          RENAME __init__.py)
else()
  install(FILES __init__.py
          DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
          COMPONENT python_bindings)
endif()
