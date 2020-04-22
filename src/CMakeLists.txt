mgis_library(MFrontGenericInterface SHARED
	  ThreadPool.cxx
	  ThreadedTaskResult.cxx
	  LibrariesManager.cxx
      Markdown.cxx
	  MatrixView.cxx
	  Variable.cxx
	  Hypothesis.cxx
	  Behaviour.cxx
	  ChangeBasis.cxx
	  State.cxx
	  BehaviourData.cxx
	  MaterialStateManager.cxx
	  MaterialDataManager.cxx
	  Integrate.cxx
	  FiniteStrainSupport.cxx)
target_include_directories(MFrontGenericInterface
   PUBLIC 
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:include>)
if(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY} Threads::Threads)
  else(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY})
  endif(Threads_FOUND)
else(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE Threads::Threads)
  endif(Threads_FOUND)
endif(UNIX)

