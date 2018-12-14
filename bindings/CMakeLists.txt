if(enable-c-bindings)
  add_subdirectory("c")
endif(enable-c-bindings)

if(enable-fortran-bindings)
  add_subdirectory("fortran")
endif(enable-fortran-bindings)

if(enable-python-bindings)
  add_subdirectory("python")
endif(enable-python-bindings)

if(enable-fenics-bindings)
  add_subdirectory("fenics")
endif(enable-fenics-bindings)
