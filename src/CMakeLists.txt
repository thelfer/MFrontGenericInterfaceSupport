mgis_library(MFrontGenericInterface SHARED
	  Raise.cxx
	  ThreadPool.cxx
	  ThreadedTaskResult.cxx
	  LibrariesManager.cxx
      Markdown.cxx
      MaterialProperty.cxx
	  RotationMatrix.cxx
	  Variable.cxx
	  Hypothesis.cxx
	  BehaviourDescription.cxx
	  Behaviour.cxx
	  State.cxx
	  BehaviourData.cxx
	  MaterialStateManager.cxx
	  MaterialDataManager.cxx
	  Integrate.cxx
	  FiniteStrainSupport.cxx
      Model.cxx)
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

