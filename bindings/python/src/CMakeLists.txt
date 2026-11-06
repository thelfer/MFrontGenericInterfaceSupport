# the mgis_python_module function is defined in
# cmake/modules/mgis-python-modules.cmake

set(MGIS_MODULE_SOURCES)
list(APPEND MGIS_MODULE_SOURCES 
 mgis-module.cxx
 ThreadPool.cxx)
if(MGIS_HAVE_TFEL)
  list(APPEND MGIS_MODULE_SOURCES
       Database.cxx)
endif(MGIS_HAVE_TFEL)
mgis_python_module(mgis _mgis ${MGIS_MODULE_SOURCES})

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

if(enable-mgis-function)
  mgis_python_module(mgis_function function
    function-module.cxx
    NumPySupport.cxx
    Function.cxx)
endif(enable-mgis-function)
