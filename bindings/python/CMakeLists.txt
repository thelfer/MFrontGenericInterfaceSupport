# store the location of the python modules in the MGIS_PYTHON_MODULES_INSTALL_DIRECTORY variable
if(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        bin/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES})
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/mgis_${MGIS_SUFFIX_FOR_PYTHON_MODULES}")
  endif(WIN32)
else(MGIS_APPEND_SUFFIX)
  if(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "bin/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/mgis")
  else(WIN32)
    set(MGIS_PYTHON_MODULES_INSTALL_DIRECTORY
        "lib${LIB_SUFFIX}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/mgis")
  endif(WIN32)
endif(MGIS_APPEND_SUFFIX)

if(enable-pybind11)
  add_subdirectory(pybind11)
else(enable-pybind11)
  add_subdirectory(boost)
endif(enable-pybind11)
