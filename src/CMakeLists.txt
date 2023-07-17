mgis_library(MFrontGenericInterface
	  Raise.cxx
	  ThreadPool.cxx
	  ThreadedTaskResult.cxx
	  LibrariesManager.cxx
      Markdown.cxx
	  RotationMatrix.cxx
	  Variable.cxx
	  Hypothesis.cxx
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
   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:include>)
if(enable-static)
  target_include_directories(MFrontGenericInterface-static
    PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
endif(enable-static)
if(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY} Threads::Threads)
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE ${MGIS_DL_LIBRARY} Threads::Threads)
    endif(enable-static)
  else(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE ${MGIS_DL_LIBRARY})
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE ${MGIS_DL_LIBRARY})
    endif(enable-static)
  endif(Threads_FOUND)
else(UNIX)
  if(Threads_FOUND)
    target_link_libraries(MFrontGenericInterface
      PRIVATE Threads::Threads)
    if(enable-static)
      target_link_libraries(MFrontGenericInterface-static
        PRIVATE Threads::Threads)
    endif(enable-static)
  endif(Threads_FOUND)
endif(UNIX)

