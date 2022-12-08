function(add_mfront_material_properties_generated_source lib file)
  set(mfront_file   "${CMAKE_CURRENT_SOURCE_DIR}/${file}.mfront")
  set(mfront_output "src/${file}-generic.cxx")
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(mfront_flags "--debug")
  else(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(mfront_flags "")
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_custom_command(
      OUTPUT  "${mfront_output}"
      COMMAND "${MFRONT}"
      ARGS    "--search-path=${PROJECT_SOURCE_DIR}/mfront/tests/properties"
      ARGS    "${mfront_flags}" "--interface=generic" "${mfront_file}"
      DEPENDS "${mfront_file}"
      COMMENT "mfront source ${mfront_file}")
  set(${lib}_SOURCES ${mfront_output} ${${lib}_SOURCES} PARENT_SCOPE)
endfunction(add_mfront_material_properties_generated_source)

function(mfront_material_properties_check_library lib)
  if(${ARGC} LESS 1)
    message(FATAL_ERROR "mfront_material_properties_library : no source specified")
  endif(${ARGC} LESS 1)
  foreach(source ${ARGN})
    add_mfront_material_properties_generated_source(${lib} ${interface} ${source})
  endforeach(source)
  add_library(${lib} SHARED EXCLUDE_FROM_ALL
    ${${lib}_SOURCES}
    ${${lib}_ADDITIONAL_SOURCES})
  target_include_directories(${lib}
      PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/include"
      PRIVATE "${TFEL_INCLUDE_PATH}")
  set_target_properties(${lib} PROPERTIES
      COMPILE_FLAGS "-DMFRONT_COMPILING")
  add_dependencies(check ${lib})
endfunction(mfront_material_properties_check_library)
