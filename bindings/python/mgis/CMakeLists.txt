add_subdirectory(fenics)

install(FILES __init__.py
        DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
        COMPONENT python_bindings)
