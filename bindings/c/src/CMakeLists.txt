mgis_library(MFrontGenericInterface-c SHARED
  Status.cxx
  Behaviour.cxx
  BehaviourData.cxx
  BehaviourDataView.cxx
  Integrate.cxx)
target_include_directories(MFrontGenericInterface-c
	PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(MFrontGenericInterface-c
	PUBLIC "${CMAKE_SOURCE_DIR}/bindings/c/include")
target_link_libraries(MFrontGenericInterface-c
   PRIVATE MFrontGenericInterface)
