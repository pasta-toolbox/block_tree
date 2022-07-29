cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(block_tree VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 ")
# custom target for libsais
# set the output destination
set(LIBSAIS ${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais/libsais.a)
set(LIBSAIS64 ${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais/libsais64.a)
# create a custom target called build_libsais that is part of ALL
# and will run each time you type make
add_custom_target(build_libsais ALL
COMMAND ${CMAKE_MAKE_PROGRAM}
WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais/
COMMENT "Original libsais makefile target")
# now create an imported static target
add_library(libsais STATIC IMPORTED)
add_library(libsais64 STATIC IMPORTED)

# Import target "libsais" for configuration ""
set_property(TARGET libsais APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_property(TARGET libsais64 APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(libsais PROPERTIES
        IMPORTED_LOCATION_NOCONFIG "${LIBSAIS}")
set_target_properties(libsais64 PROPERTIES
        IMPORTED_LOCATION_NOCONFIG "${LIBSAIS64}")
# now you can use libsaus as if it were a regular cmake built target in your project
add_dependencies(libsais LIBSAIS)
add_dependencies(libsais64 LIBSAIS64)
add_executable(foo foo.cpp)
add_subdirectory(includes)
target_link_libraries(foo libsais)
target_link_libraries(foo libsais64)
target_link_libraries(
        foo
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/malloc_count/malloc_count.o
)
add_subdirectory(extlib/tlx)
add_subdirectory(extlib/bit_vector)
target_link_libraries(foo pasta_bit_vector)
target_link_libraries(foo tlx)
target_include_directories(foo PUBLIC extlib/bit_vector/include/pasta/bit_vector)
target_include_directories(foo PUBLIC extlib/tlx/tlx)
target_include_directories(foo PUBLIC extlib/bit_vector/include/pasta/bit_vector/support)
target_include_directories(foo PUBLIC includes/utils)
target_include_directories(foo PUBLIC includes/bv_block_tree)
target_include_directories(foo PUBLIC extlib/libsais/src)
target_include_directories(foo PUBLIC extlib/malloc_count)