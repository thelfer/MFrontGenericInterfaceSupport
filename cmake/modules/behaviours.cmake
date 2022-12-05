macro(add_mfront_dependency name file)
  set(mfront_file    "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront")
  set(mfront_output1 "${CMAKE_CURRENT_BINARY_DIR}/src/${file}-mfront.cxx")
  set(mfront_output2 "${CMAKE_CURRENT_BINARY_DIR}/include/${file}-mfront.hxx")
  add_custom_command(
    OUTPUT  "${mfront_output1}"
    OUTPUT  "${mfront_output2}"
    COMMAND "${MFRONT}"
    ARGS    "--search-path=${CMAKE_CURRENT_SOURCE_DIR}"
    ARGS    "--interface=mfront" "${mfront_file}"
    DEPENDS "${mfront_file}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "mfront source ${mfront_file}")
  set(${name}_dependencies_SOURCES
      ${mfront_output1} ${mfront_output2}
      ${${name}_dependencies_SOURCES})
endmacro(add_mfront_dependency)

function(add_mfront_behaviour_sources lib  file)
  set(mfront_file   "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront")
  set(mfront_output1 "src/${file}.cxx")
  set(mfront_output2 "src/${file}-generic.cxx")
  add_custom_command(
    OUTPUT  "${mfront_output1}"
    OUTPUT  "${mfront_output2}"
    COMMAND "${MFRONT}"
    ARGS    "--interface=generic"
    ARGS     "${mfront_file}"
    DEPENDS "${mfront_file}"
    COMMENT "mfront source ${mfront_file}")
  set(${lib}_SOURCES ${mfront_output1} ${mfront_output2}
    ${${lib}_SOURCES} PARENT_SCOPE)
endfunction(add_mfront_behaviour_sources)

function(mfront_behaviours_check_library name)
  unset(SOURCES)    # Potentially set on the caller side.
  set ( _CMD SOURCES )
  set ( _SOURCES )
  set ( _DEPENDENCIES )
  foreach ( _ARG ${ARGN})
    if ( ${_ARG} MATCHES SOURCES )
      set ( _CMD SOURCES )
    elseif ( ${_ARG} MATCHES DEPENDENCIES )
      set ( _CMD DEPENDENCIES )
    else ()
      if ( ${_CMD} MATCHES SOURCES )
        list ( APPEND _SOURCES "${_ARG}" )
      elseif ( ${_CMD} MATCHES DEPENDENCIES )
        list ( APPEND _DEPENDENCIES "${_ARG}" )
      endif ()
    endif ()
  endforeach ()
  list(LENGTH _SOURCES _SOURCES_LENGTH )
  if(${_SOURCES_LENGTH} LESS 1)
    message(FATAL_ERROR "mfront_behaviours_check_library : no source specified")
  endif(${_SOURCES_LENGTH} LESS 1)
  # treating dependencies
  foreach(dep ${_DEPENDENCIES})
    add_mfront_dependency(${name} ${dep})
  endforeach(dep ${_DEPENDENCIES})
  include_directories("${TFEL_INCLUDE_PATH}")
  foreach(source ${_SOURCES})
    add_mfront_behaviour_sources(${name} ${source})
  endforeach(source ${_SOURCES})
  if(${name}_dependencies_SOURCES)
    message("dependencies : ${${name}_dependencies_SOURCES}")
  endif(${name}_dependencies_SOURCES)
  foreach(deps ${${name}_dependencies_SOURCES})
    set(${name}_SOURCES ${deps} ${${name}_SOURCES})
  endforeach(deps ${${name}_dependencies_SOURCES})
  list(LENGTH ${name}_SOURCES nb_sources)
  if(nb_sources GREATER 0)
    message(STATUS "Adding library : ${name} (${${name}_SOURCES})")
    add_library(${name} SHARED EXCLUDE_FROM_ALL ${${name}_SOURCES})
    add_dependencies(check ${name})
    target_include_directories(${name}
      PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/${interface}/include"
      PRIVATE "${TFEL_INCLUDE_PATH}")
    if(WIN32)
      if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
	set_target_properties(${name}
	  PROPERTIES LINK_FLAGS "-Wl,--kill-at -Wl,--no-undefined")
      endif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    endif(WIN32)
    target_link_libraries(${name} PUBLIC ${TFEL_MFRONT_LIBRARIES})
  else(nb_sources GREATER 0)
    message(STATUS "No sources selected for library ${name}")
  endif(nb_sources GREATER 0)
endfunction(mfront_behaviours_check_library)
