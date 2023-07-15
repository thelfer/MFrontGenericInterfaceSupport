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
    # set_target_properties(${name}-static PROPERTIES OUTPUT_NAME "${name}-static")
    # set_target_properties(${name}-static PROPERTIES COMPILE_FLAGS "-D${name}_EXPORTS -DMGIS_STATIC_BUILD")
    if(WIN32)
      install(TARGETS ${name}-static EXPORT ${name}-static DESTINATION bin)
    else(WIN32)
      install(TARGETS ${name}-static EXPORT ${name}-static DESTINATION lib${LIB_SUFFIX})
    endif(WIN32)
    # install(EXPORT ${name}-static DESTINATION ${export_install_path}
    #         NAMESPACE mgis:: FILE ${name}Config.cmake)
  endif(enable-static)
endfunction(mgis_library)
