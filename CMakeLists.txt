cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(block_tree LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 ")

add_executable(foo foo.cpp includes/utils/print_util.h)
add_subdirectory(includes)
target_link_libraries(
        foo
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/malloc_count/malloc_count.o
)
add_subdirectory(extlib/libsais_cmake)
add_subdirectory(extlib/tlx)
add_subdirectory(extlib/bit_vector)
target_link_libraries(foo pasta_bit_vector)
target_link_libraries(foo tlx)
target_link_libraries(foo sais)
target_link_libraries(foo sais64)
target_include_directories(foo PUBLIC extlib/bit_vector/include/pasta/bit_vector)
target_include_directories(foo PUBLIC extlib/tlx/tlx)
target_include_directories(foo PUBLIC extlib/bit_vector/include/pasta/bit_vector/support)
target_include_directories(foo PUBLIC includes/utils)
target_include_directories(foo PUBLIC includes/bv_block_tree)
target_include_directories(foo PUBLIC extlib/libsais_cmake/src)
target_include_directories(foo PUBLIC extlib/malloc_count)
