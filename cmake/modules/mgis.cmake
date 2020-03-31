function(mgis_header dir file)
  install(FILES ${dir}/${file}
    DESTINATION "include/${dir}")
endfunction(mgis_header)

function(mgis_library name)
  if(${ARGC} LESS 2)
    message(FATAL_ERROR "mgis_library_internal : no source specified")
  endif(${ARGC} LESS 2)
  add_library(${name} SHARED ${ARGN})
  if(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION bin)
  else(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION lib${LIB_SUFFIX})
  endif(WIN32)
  if(MGIS_APPEND_SUFFIX)
    set(export_install_path "share/mgis-${MGIS_SUFFIX}/cmake")
  else(MGIS_APPEND_SUFFIX)
    set(export_install_path "share/mgis/cmake")
  endif(MGIS_APPEND_SUFFIX)
  install(EXPORT ${name} DESTINATION ${export_install_path}
          NAMESPACE mgis:: FILE ${name}Config.cmake)
  if(enable-static)
    add_library(${name}-static STATIC ${ARGN})
    set_target_properties(${name}-static PROPERTIES OUTPUT_NAME "${name}")
    # Now the library target "${name}-static" will be named "${name}.lib"
    # with MS tools.
    # This conflicts with the "${name}.lib" import library corresponding
    # to "${name}.dll",
    # so we add a "lib" prefix (which is default on other platforms
    # anyway):
    set_target_properties(${name}-static PROPERTIES PREFIX "lib")
    # Help CMake 2.6.x and lower (not necessary for 2.8 and above, but
    # doesn't hurt):
    set_target_properties(${name}        PROPERTIES CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(${name}-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
    set_target_properties(${name}-static PROPERTIES COMPILE_FLAGS "-D${name}_EXPORTS -DMGIS_STATIC_BUILD")
    if(WIN32)
      install(TARGETS ${name}-static DESTINATION bin)
    else(WIN32)
      install(TARGETS ${name}-static DESTINATION lib${LIB_SUFFIX})
    endif(WIN32)
    install(EXPORT ${name}-static DESTINATION ${export_install_path}
            NAMESPACE mgis:: FILE ${name}Config.cmake)
  endif(enable-static)
endfunction(mgis_library)
