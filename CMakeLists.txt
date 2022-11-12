cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(block_tree LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp -lm -lstdc++")
# 
#-fsanitize=address
#add_executable(foo benchmark/foo.cpp includes/utils/print_util.h)
add_library(block_tree INTERFACE)
# Options when compiling pasta::bit_vector
## Build tests
option(BLOCK_TREE_BUILD_TESTS
        "Build blocktree's tests." OFF)
option(BLOCK_TREE_BUILD_BENCHMARKS
        "Build blocktree's benchmarks." OFF)
target_include_directories(block_tree INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
# Optional test
if(BLOCK_TREE_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.12.1
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    add_subdirectory(tests)
    include(GoogleTest)
endif()
#if(BLOCK_TREE_BUILD_BENCHMARKS)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmark)
#endif()

add_subdirectory(includes)
# add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/malloc_count_cmake)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais_cmake)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/bit_vector)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/tlx)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/sdsl)
target_link_libraries(block_tree INTERFACE pasta_bit_vector)
target_link_libraries(block_tree INTERFACE sdsl)
target_link_libraries(block_tree INTERFACE tlx)
target_link_libraries(block_tree INTERFACE sais)
target_link_libraries(block_tree INTERFACE sais64)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/bit_vector/include/pasta/bit_vector)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/tlx/tlx)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/bit_vector/include/pasta/bit_vector/support)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/sdsl/include/sdsl)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/includes/utils)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/includes/bv_block_tree)
target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais_cmake/src)
# target_include_directories(block_tree INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/extlib/malloc_count_cmake)

