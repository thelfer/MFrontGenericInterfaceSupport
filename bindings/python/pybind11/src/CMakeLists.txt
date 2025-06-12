function(mgis_python_module fullname name)
  if(${ARGC} LESS 1)
    message(FATAL_ERROR "mgis_python_module: no source specified")
  endif(${ARGC} LESS 1)
  add_library(py_${fullname} MODULE ${ARGN})
  target_include_directories(py_${fullname} PRIVATE
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/bindings/python/pybind11/include>)
  if(TARGET pybind11::lto)
    target_link_libraries(py_${fullname} PRIVATE
      MFrontGenericInterface
      pybind11::module pybind11::lto pybind11::windows_extras)
  else(TARGET pybind11::lto)
    target_link_libraries(py_${fullname} PRIVATE
      MFrontGenericInterface
      pybind11::module pybind11::windows_extras)
  endif(TARGET pybind11::lto)
  pybind11_extension(py_${fullname})
  if(NOT MSVC AND NOT ${CMAKE_BUILD_TYPE} MATCHES Debug|RelWithDebInfo)
    # Strip unnecessary sections of the binary on Linux/macOS
    pybind11_strip(py_${fullname})
  endif()
  set_target_properties(py_${fullname} PROPERTIES PREFIX "")
  set_target_properties(py_${fullname} PROPERTIES OUTPUT_NAME ${name})
  set_target_properties(py_${fullname} PROPERTIES
    CXX_VISIBILITY_PRESET "hidden"
    CUDA_VISIBILITY_PRESET "hidden")
  # installation
  install(TARGETS py_${fullname}
    DESTINATION ${MGIS_PYTHON_MODULES_INSTALL_DIRECTORY}
    COMPONENT python_bindings)
endfunction(mgis_python_module)

mgis_python_module(mgis _mgis
 mgis-module.cxx
 ThreadPool.cxx)

mgis_python_module(mgis_material_property material_property
 material_property-module.cxx
 MaterialProperty.cxx)

mgis_python_module(mgis_behaviour behaviour
 behaviour-module.cxx
 NumPySupport.cxx
 Hypothesis.cxx
 Variable.cxx
 BehaviourDescription.cxx
 Behaviour.cxx
 State.cxx
 BehaviourData.cxx
 BehaviourDataView.cxx
 MaterialDataManager.cxx
 MaterialStateManager.cxx
 Integrate.cxx
 FiniteStrainSupport.cxx)

mgis_python_module(mgis_model model
 model-module.cxx
 Model.cxx)


