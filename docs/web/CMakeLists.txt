if(MGIS_HAVE_PANDOC)
  function(mgis_pandoc_generate_html_page file)
    ADD_CUSTOM_COMMAND(
      OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/${file}.md
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/mgis-template.html
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/css/main.css
      COMMAND   ${PANDOC}
      ARGS      -f markdown-markdown_in_html_blocks+tex_math_single_backslash
      ARGS      --filter pandoc-crossref -M crossrefYaml=${CMAKE_CURRENT_SOURCE_DIR}/pandoc-crossref.yaml
      ARGS      --bibliography=${CMAKE_CURRENT_SOURCE_DIR}/bibliography.bib --filter pandoc-citeproc
      ARGS      --template=${CMAKE_CURRENT_SOURCE_DIR}/mgis-template.html
      ARGS      --mathjax
      ARGS      --email-obfuscation=javascript
      ARGS       --default-image-extension=svg
      ARGS      ${ARGN}
      ARGS      ${CMAKE_CURRENT_SOURCE_DIR}/${file}.md -o ${file}.html)
    add_custom_target(${file}-html ALL DEPENDS ${file}.html)
    add_dependencies(website ${file}-html)
    if(MGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
        COMPONENT website)
    else(MGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mgis/web
        COMPONENT website)
    endif(MGIS_APPEND_SUFFIX)
  endfunction(mgis_pandoc_generate_html_page)
endif(MGIS_HAVE_PANDOC)

if(MGIS_HAVE_JUPYTER_NBCONVERT)
  function(mgis_jupyter_nbconvert_generate_html_page file)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${file}.ipynb
                   ${CMAKE_CURRENT_BINARY_DIR}/${file}.ipynb COPYONLY)
    ADD_CUSTOM_COMMAND(
      OUTPUT    ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
      DEPENDS   ${CMAKE_CURRENT_BINARY_DIR}/${file}.ipynb
      COMMAND   ${JUPYTER_NBCONVERT}
      ARGS      --to html
      ARGS      ${CMAKE_CURRENT_BINARY_DIR}/${file}.ipynb)
    add_custom_target(${file}-html ALL DEPENDS ${file}.html)
    add_dependencies(website ${file}-html)
    if(MGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
        COMPONENT website)
    else(MGIS_APPEND_SUFFIX)
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${file}.html
        DESTINATION share/doc/mgis/web
        COMPONENT website)
    endif(MGIS_APPEND_SUFFIX)
  endfunction(mgis_jupyter_nbconvert_generate_html_page)
endif(MGIS_HAVE_JUPYTER_NBCONVERT)

add_custom_target(website)
add_dependencies(doc website)

if(MGIS_HAVE_PANDOC)
  mgis_pandoc_generate_html_page(index)
  mgis_pandoc_generate_html_page(install)
  mgis_pandoc_generate_html_page(bindings-cxx "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(FEniCSBindings)
  mgis_pandoc_generate_html_page(release-notes-1.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.2 "--toc" "--toc-depth=3")
endif(MGIS_HAVE_PANDOC)

if(MGIS_HAVE_JUPYTER_NBCONVERT)
  mgis_jupyter_nbconvert_generate_html_page(mgis_fenics_multiphase_model)
endif(MGIS_HAVE_JUPYTER_NBCONVERT)


if(MGIS_APPEND_SUFFIX)
  install(DIRECTORY css
    DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
#  install(DIRECTORY js css doc img media downloads
#    DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
#    COMPONENT website
#    PATTERN ".svn" EXCLUDE)
#  install(FILES 404.html apple-touch-icon-precomposed.png 
#    CHANGELOG.md CONTRIBUTING.md crossdomain.xml favicon.ico
#    humans.txt LICENSE.md README.md robots.txt TODO.md
#    DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
#    COMPONENT website)
else(MGIS_APPEND_SUFFIX)
  install(DIRECTORY css
    DESTINATION share/doc/mgis/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
#  install(DIRECTORY js css doc img media downloads
#    DESTINATION share/doc/mgis/web
#    COMPONENT website
#    PATTERN ".svn" EXCLUDE)
#  install(FILES 404.html apple-touch-icon-precomposed.png 
#    CHANGELOG.md CONTRIBUTING.md crossdomain.xml favicon.ico
#    humans.txt LICENSE.md README.md robots.txt TODO.md
#    DESTINATION share/doc/mgis/web
#    COMPONENT website)
endif(MGIS_APPEND_SUFFIX)
