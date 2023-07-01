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

# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20  -g -fopenmp -lm -lstdc++")
# 
#-g -fsanitize=address
#add_executable(foo benchmark/foo.cpp includes/utils/print_util.h)
# Options when compiling pasta::bit_vector
## Build tests
option(PASTA_BLOCK_TREE_BUILD_TESTS
        "Build blocktree's tests." OFF)
option(PASTA_BLOCK_TREE_BUILD_BENCHMARKS
  "Build blocktree's benchmarks." OFF)

# Optional test
if(PASTA_BLOCK_TREE_BUILD_TESTS)
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
if(PASTA_BLOCK_TREE_BUILD_BENCHMARKS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/benchmark)
endif()

set(LIBSAIS_USE_OPENMP ON CACHE BOOL "Use OpenMP for parallelization of libsais" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/libsais)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/bit_vector)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extlib/tlx)

add_library(pasta_block_tree INTERFACE)
target_include_directories(pasta_block_tree INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(pasta_block_tree INTERFACE
  libsais
  pasta_bit_vector
  sdsl
  tlx)

################################################################################
