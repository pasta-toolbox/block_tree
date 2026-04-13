################################################################################
# This file is part of pasta::block_tree
#
# Copyright (C) 2022 Daniel Meyer
# Copyright (C) 2023 Florian Kurpicz <florian@kurpicz.org>
#
# pasta::block_tree is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pasta::block_tree is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pasta::block_tree.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(pasta_block_tree)

# Generate compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

## Build tests
option(PASTA_BLOCK_TREE_BUILD_TESTS
        "Build blocktree's tests." OFF)
option(PASTA_BLOCK_TREE_BUILD_EXAMPLES
        "Build blocktree's benchmarks." OFF)
option(PASTA_BLOCK_TREE_DEBUG
        "Print debug information" OFF)
option(PASTA_BLOCK_TREE_BENCH
        "Enable outputting benchmark information" OFF)

include(ExternalProject)
include(FetchContent)

FetchContent_Declare(
  tlx
  GIT_REPOSITORY https://github.com/tlx/tlx.git
  GIT_TAG        b6af589 #release 0.6.1
)

FetchContent_Declare(
  pasta_bit_vector
  GIT_REPOSITORY https://github.com/pasta-toolbox/bit_vector.git
  GIT_TAG        b4798d5 #main
)

FetchContent_MakeAvailable(tlx pasta_bit_vector)

# Optional test
if (PASTA_BLOCK_TREE_BUILD_TESTS)
    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG release-1.12.1
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(tests)
endif ()
if (PASTA_BLOCK_TREE_BUILD_EXAMPLES)
    add_executable(block_tree_construction
            examples/block_tree_construction.cpp)
    target_link_libraries(block_tree_construction
            pasta_block_tree)

    add_executable(build_bt
            examples/build_bt.cpp)
    target_link_libraries(build_bt
            pasta_block_tree)
    if (PASTA_BLOCK_TREE_DEBUG)
        target_compile_definitions(build_bt PRIVATE BT_INSTRUMENT)
        target_compile_definitions(build_bt PRIVATE BT_DBG)
    endif ()

    if (PASTA_BLOCK_TREE_BENCH)
        target_compile_definitions(build_bt PRIVATE BT_INSTRUMENT)
        target_compile_definitions(build_bt PRIVATE BT_BENCH)
    endif ()
endif ()

set(LIBSAIS_USE_OPENMP ON CACHE BOOL "Use OpenMP for parallelization of libsais" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais)

set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/sdsl-lite/external/libdivsufsort)

add_library(waitfree-mpsc-queue
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/waitfree-mpsc-queue/mpsc.c)
target_include_directories(waitfree-mpsc-queue PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/waitfree-mpsc-queue)

add_library(jiffy INTERFACE)
target_include_directories(jiffy INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/Jiffy)
add_library(jiffy1 INTERFACE)
target_include_directories(jiffy1 INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/Jiffy-1)

add_library(pasta_block_tree INTERFACE)
target_include_directories(pasta_block_tree INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include)

# sdsl
file(GLOB sdsl_sources ${CMAKE_CURRENT_SOURCE_DIR}/extlib/sdsl-lite/lib/*.cpp)
add_library(sdsl STATIC ${sdsl_sources})
target_include_directories(sdsl SYSTEM PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/sdsl-lite/include)
target_include_directories(sdsl SYSTEM PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}/extlib/sdsl-lite/external/libdivsufsort/include)
target_link_libraries(sdsl PUBLIC divsufsort64)
set_target_properties(sdsl PROPERTIES COMPILE_FLAGS "-w")

target_link_libraries(pasta_block_tree INTERFACE
        libsais
        pasta_bit_vector
        tlx
        waitfree-mpsc-queue
        sdsl
        #jiffy
        jiffy1
        #        libzstd_static
)

target_include_directories(pasta_block_tree INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/growt)
target_include_directories(pasta_block_tree INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/parallel-hashmap/parallel_hashmap)
target_include_directories(pasta_block_tree SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/parlayhash/include)
target_include_directories(pasta_block_tree SYSTEM INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/extlib/unordered_dense/include)

################################################################################
