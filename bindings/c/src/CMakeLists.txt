mgis_library(MFrontGenericInterface-c SHARED
  Status.cxx
  ThreadPool.cxx
  Behaviour.cxx
  State.cxx
  BehaviourData.cxx
  BehaviourDataView.cxx
  MaterialStateManager.cxx
  MaterialDataManager.cxx
  Integrate.cxx)
target_include_directories(MFrontGenericInterface-c
  PUBLIC "${CMAKE_SOURCE_DIR}/include"
  PUBLIC "${CMAKE_SOURCE_DIR}/bindings/c/include")
target_link_libraries(MFrontGenericInterface-c
   PRIVATE MFrontGenericInterface)
