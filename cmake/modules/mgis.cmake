include(CMakePackageConfigHelpers)

if(MGIS_APPEND_SUFFIX)
  set(mgis_export_install_path "share/mgis-${MGIS_SUFFIX}/cmake")
else(MGIS_APPEND_SUFFIX)
  set(mgis_export_install_path "share/mgis/cmake")
endif(MGIS_APPEND_SUFFIX)

function(mgis_header dir file)
  install(FILES ${dir}/${file}
    DESTINATION "include/${dir}")
endfunction(mgis_header)

function(mgis_library name)
  if(${ARGC} LESS 2)
    message(FATAL_ERROR "mgis_library_internal : no source specified")
  endif(${ARGC} LESS 2)
  add_library(${name} SHARED ${ARGN})
  add_library(mgis::${name} ALIAS ${name})
  if(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION bin)
  else(WIN32)
    install(TARGETS ${name} EXPORT ${name}
            DESTINATION lib${LIB_SUFFIX})
  endif(WIN32)
  install(EXPORT ${name} DESTINATION ${mgis_export_install_path}
          EXPORT_LINK_INTERFACE_LIBRARIES
          NAMESPACE mgis:: FILE ${name}Targets.cmake)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name}Config.cmake.in)
    set(_package_config_file ${CMAKE_CURRENT_SOURCE_DIR}/${name}Config.cmake.in)
  else()
    set(_package_config_file ${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake.in)
    file(WRITE ${_package_config_file}
         "@PACKAGE_INIT@\n"
         "\n"
         "include(\"\${CMAKE_CURRENT_LIST_DIR}/${name}Targets.cmake\")")
  endif()
  # generate the config file that includes the exports
  configure_package_config_file(${_package_config_file}
                                "${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake"
                                INSTALL_DESTINATION ${mgis_export_install_path}
                                NO_SET_AND_CHECK_MACRO
                                NO_CHECK_REQUIRED_COMPONENTS_MACRO)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${name}Config.cmake"
          DESTINATION ${mgis_export_install_path})
  if(enable-static)
    add_library(${name}-static STATIC ${ARGN})
    # set_target_properties(${name}-static PROPERTIES OUTPUT_NAME "${name}-static")
    # set_target_properties(${name}-static PROPERTIES COMPILE_FLAGS "-D${name}_EXPORTS -DMGIS_STATIC_BUILD")
    if(WIN32)
      install(TARGETS ${name}-static EXPORT ${name}-static DESTINATION bin)
    else(WIN32)
      install(TARGETS ${name}-static EXPORT ${name}-static DESTINATION lib${LIB_SUFFIX})
    endif(WIN32)
    # install(EXPORT ${name}-static DESTINATION ${mgis_export_install_path}
    #         NAMESPACE mgis:: FILE ${name}Config.cmake)
  endif(enable-static)
endfunction(mgis_library)
