mgis_library(MFrontGenericInterface SHARED
  LibrariesManager.cxx
  Variable.cxx
  Hypothesis.cxx
  Description.cxx)
target_include_directories(MFrontGenericInterface
	PUBLIC "${CMAKE_SOURCE_DIR}/include")
if(UNIX)
  target_link_libraries(MFrontGenericInterface
    PRIVATE ${MFRONT_DL_LIBRARY})
endif(UNIX)

