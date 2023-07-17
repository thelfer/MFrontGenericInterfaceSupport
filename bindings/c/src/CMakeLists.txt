mgis_library(MFrontGenericInterface-c
  Status.cxx
  ThreadPool.cxx
  Hypothesis.cxx
  Variable.cxx
  Behaviour.cxx
  State.cxx
  BehaviourData.cxx
  BehaviourDataView.cxx
  MaterialStateManager.cxx
  MaterialDataManager.cxx
  Integrate.cxx
  Model.cxx)
target_include_directories(MFrontGenericInterface-c
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/bindings/c/include>
  PUBLIC $<INSTALL_INTERFACE:include>)
target_link_libraries(MFrontGenericInterface-c
   PRIVATE MFrontGenericInterface)
if(enable-static)
  target_include_directories(MFrontGenericInterface-c-static
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/bindings/c/include>
    PUBLIC $<INSTALL_INTERFACE:include>)
  target_link_libraries(MFrontGenericInterface-c-static
     PRIVATE MFrontGenericInterface)
endif(enable-static)
