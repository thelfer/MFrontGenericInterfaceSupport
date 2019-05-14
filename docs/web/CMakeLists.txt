macro(mgis_generate_html_page file)
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
endmacro(mgis_generate_html_page)

add_custom_target(website)
add_dependencies(doc website)

mgis_generate_html_page(index)
mgis_generate_html_page(bindings-cxx "--toc" "--toc-depth=3")
mgis_generate_html_page(FEniCSBindings)

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
