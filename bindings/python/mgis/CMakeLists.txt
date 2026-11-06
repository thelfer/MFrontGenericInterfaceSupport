add_subdirectory(fenics)

if(MGIS_HAVE_TFEL AND TFEL_PYTHON_BINDINGS)
  set(mgis_python_init_script __init__-have-tfel.py)
else()
  set(mgis_python_init_script __init__.py)
endif()

install(FILES ${mgis_python_init_script}
        DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
        COMPONENT python_bindings
        RENAME __init__.py)
# PEP 561 marker: the package installs the .pyi stubs of its pybind11
# modules, which type checkers only use if this marker is present
install(FILES py.typed
        DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
        COMPONENT python_bindings)
# In the build tree, the package directory where the modules are placed
# needs the __init__.py script so that the package can be imported, which
# pybind11-stubgen requires to generate the stubs.
configure_file(${mgis_python_init_script}
               ${MGIS_PYTHON_MODULES_BUILD_DIR}/${MGIS_PYTHON_PACKAGE_NAME}/__init__.py
               COPYONLY)
