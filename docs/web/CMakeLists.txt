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
      ARGS      --bibliography=${CMAKE_CURRENT_SOURCE_DIR}/bibliography.bib --citeproc
      ARGS      --template=${CMAKE_CURRENT_SOURCE_DIR}/mgis-template.html --css=${CMAKE_CURRENT_SOURCE_DIR}/css/main.css
      ARGS      --mathjax
      ARGS      --highlight-style=tango
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
if(enable-website)
add_dependencies(doc website)
endif(enable-website)
  
if(MGIS_HAVE_PANDOC)
  mgis_pandoc_generate_html_page(index "-c css/slideshow.css")
  mgis_pandoc_generate_html_page(install)
  mgis_pandoc_generate_html_page(contributing)
  mgis_pandoc_generate_html_page(bindings-cxx "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(bindings-python-small_strain "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(bindings-python-finite_strain "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.0.2 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.1.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.2 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.2.1 "--number-sections" "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.2.2 "--number-sections" "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-1.2.3 "--number-sections" "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-2.0 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-2.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-2.2 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-2.2.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-3.0 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-3.0.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(release-notes-3.1 "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(orthotropic-behaviours "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(FEniCSBindings)
  mgis_pandoc_generate_html_page(mgis_fenics "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_nonlinear_heat_transfer "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_nonlinear_heat_transfer_3D "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_small_strain_elastoplasticity "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_finite_strain_elastoplasticity "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_multiphase_model "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_phase_field "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_heat_equation_phase_change "--toc" "--toc-depth=3")
  mgis_pandoc_generate_html_page(mgis_fenics_monolithic_transient_thermoelasticity "--toc" "--toc-depth=3")
endif(MGIS_HAVE_PANDOC)

if(MGIS_HAVE_JUPYTER_NBCONVERT)
#  mgis_jupyter_nbconvert_generate_html_page(mgis_fenics_multiphase_model)
endif(MGIS_HAVE_JUPYTER_NBCONVERT)

if(MGIS_APPEND_SUFFIX)
  install(DIRECTORY css img js mfront ipynb
    DESTINATION share/doc/mgis-${MGIS_SUFFIX}/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
else(MGIS_APPEND_SUFFIX)
  install(DIRECTORY css img js mfront ipynb
    DESTINATION share/doc/mgis/web
    COMPONENT website
    PATTERN ".svn" EXCLUDE)
endif(MGIS_APPEND_SUFFIX)
