# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(linear_cell_complex_experimentations)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
find_package( TBB REQUIRED )
include(CGAL_TBB_support)

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math ckeck.
# add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

create_single_source_cgal_program("iterate_over_vertex.cpp")
create_single_source_cgal_program("multithread_dart_iteration.cpp")
create_single_source_cgal_program("concurrent_compact_container_test.cpp")
create_single_source_cgal_program("benchmark.cpp")

if(CGAL_Qt6_FOUND)
  target_link_libraries(iterate_over_vertex PUBLIC CGAL::CGAL_Basic_viewer)
endif()

target_link_libraries(concurrent_compact_container_test PUBLIC CGAL::TBB_support)
target_link_libraries(benchmark PUBLIC CGAL::TBB_support)
target_link_libraries(multithread_dart_iteration PUBLIC CGAL::TBB_support)
