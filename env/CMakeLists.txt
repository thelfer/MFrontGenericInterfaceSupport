if(UNIX)
  if(enable-python-bindings)
    configure_file(env-python.sh.in env.sh @ONLY)
    configure_file(env-python.fish.in env.fish @ONLY)
  else(enable-python-bindings) 
   configure_file(env.sh.in env.sh @ONLY)
   configure_file(env.fish.in env.fish @ONLY)
  endif(enable-python-bindings)
  if(MGIS_SUFFIX)
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/env.sh
            ${CMAKE_CURRENT_BINARY_DIR}/env.fish
            DESTINATION share/mgis-${MGIS_SUFFIX}/env
            COMPONENT core)
  else(MGIS_SUFFIX)
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/env.sh
            ${CMAKE_CURRENT_BINARY_DIR}/env.fish
            DESTINATION share/mgis/env
            COMPONENT core)
  endif(MGIS_SUFFIX)
endif(UNIX)

if(WIN32)
  if(enable-python-bindings)
    cmake_path(NATIVE_PATH MGIS_PYTHON_MODULES_INSTALL_DIRECTORY NATIVE_MGIS_PYTHON_MODULES_INSTALL_DIRECTORY)
    configure_file(env-python.ps1.in env.ps1 @ONLY)
    configure_file(env-python.bat.in env.bat @ONLY)
  else(enable-python-bindings) 
   configure_file(env.ps1.in env.ps1 @ONLY)
   configure_file(env.bat.in env.bat @ONLY)
  endif(enable-python-bindings)
  if(MGIS_SUFFIX)
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/env.ps1
            ${CMAKE_CURRENT_BINARY_DIR}/env.bat
            DESTINATION share/mgis-${MGIS_SUFFIX}/env
            COMPONENT core)
  else(MGIS_SUFFIX)
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/env.ps1
            ${CMAKE_CURRENT_BINARY_DIR}/env.bat
            DESTINATION share/mgis/env
            COMPONENT core)
  endif(MGIS_SUFFIX)
endif(WIN32)
