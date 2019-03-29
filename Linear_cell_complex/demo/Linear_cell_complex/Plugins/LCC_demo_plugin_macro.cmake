include(AddFileDependencies)
include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)

  macro(LCC_demo_plugin plugin_name plugin_implementation_base_name)
    list_split(option ARGN_TAIL ${ARGN} )
    if(NOT ${option} STREQUAL "EXCLUDE_FROM_ALL")
      if(NOT ${option} STREQUAL "NO_MOC")
        set(other_sources ${ARGN})
        set(option "")
      else()
        set(other_sources ${ARGN_TAIL})
      endif()
    else()
      set(other_sources ${ARGN_TAIL})
    endif()
    if("${option}" STREQUAL "NO_MOC") 
      set(option "")
      set(moc_file_name "")
    else()
      set(moc_file_name ${plugin_implementation_base_name}.moc)
      qt5_generate_moc("${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp"
                        "${CMAKE_CURRENT_BINARY_DIR}/${moc_file_name}")      
      add_file_dependencies(${moc_file_name}
                            "${CMAKE_CURRENT_SOURCE_DIR}/${plugin_implementation_base_name}.cpp")
    endif()

    add_library(${plugin_name} MODULE ${option} ${moc_file_name}
                ${plugin_implementation_base_name}.cpp ${other_sources})
    cgal_add_compilation_test(${plugin_name})

    add_to_cached_list( CGAL_EXECUTABLE_TARGETS ${plugin_name} )
    target_link_libraries( ${plugin_name} PUBLIC ${QT_LIBRARIES} )
    target_link_libraries( ${plugin_name} PUBLIC ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
  endmacro(LCC_demo_plugin)
