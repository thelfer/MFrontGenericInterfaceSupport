add_subdirectory(src)
add_subdirectory(modules)
if(MGIS_HAVE_MFRONT_SUPPORT)
  add_subdirectory(tests)
endif(MGIS_HAVE_MFRONT_SUPPORT)

