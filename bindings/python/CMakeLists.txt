# store the location of the python modules in the MGIS_PYTHON_MODULES_INSTALL_DIRECTORY variable
if(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        bin/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES})
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES}")
  endif(WIN32)
else(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "bin/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/mgis")
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/mgis")
  endif(WIN32)
endif(MGIS_APPEND_SUFFIX)

add_subdirectory(src)
add_subdirectory(mgis)
if(MGIS_HAVE_TFEL)
  add_subdirectory(tests)
endif(MGIS_HAVE_TFEL)
