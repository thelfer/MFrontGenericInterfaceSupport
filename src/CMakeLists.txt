mgis_library(MFrontGenericInterface SHARED
	  ThreadPool.cxx
	  ThreadedTaskResult.cxx
	  LibrariesManager.cxx
	  Variable.cxx
	  Hypothesis.cxx
	  Behaviour.cxx
	  State.cxx
	  BehaviourData.cxx
	  MaterialStateManager.cxx
	  MaterialDataManager.cxx
	  Integrate.cxx)
target_include_directories(MFrontGenericInterface
	PUBLIC "${CMAKE_SOURCE_DIR}/include")
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

