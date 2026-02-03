if(MGIS_APPEND_SUFFIX)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION share/doc/mgis-${MGIS_SUFFIX}
    COMPONENT doxygen)
  install(FILES
    ${PROJECT_SOURCE_DIR}/docs/doxygen/img/background.svg
    DESTINATION share/doc/mgis-${MGIS_SUFFIX}/html/img/
    COMPONENT doxygen)
else(MGIS_APPEND_SUFFIX)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION share/doc/mgis
    COMPONENT doxygen)
  install(FILES
    ${PROJECT_SOURCE_DIR}/docs/doxygen/img/background.svg
    DESTINATION share/doc/mgis/html/img/
    COMPONENT doxygen)
endif(MGIS_APPEND_SUFFIX)
