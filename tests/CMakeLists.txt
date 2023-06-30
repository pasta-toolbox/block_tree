################################################################################
# This file is part of pasta::block_tree
#
# Copyright (C) 2022 Daniel Meyer
#
# PaStA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PaStA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PaStA.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

macro(pasta_block_tree_build_test TESTNAME)
    string(REPLACE "/" "_" TESTNAME_REPLACED "${TESTNAME}")
    add_executable(${TESTNAME_REPLACED}
            ${TESTNAME}.cpp)
    target_link_libraries(${TESTNAME_REPLACED} PRIVATE
            pasta_block_tree
            GTest::Main)
    include_directories(${TESTNAME_REPLACED} PRIVATE ${gtest_SOURCE_DIR}/include)
    add_test(
            NAME ${TESTNAME_REPLACED}
            COMMAND ${TESTNAME_REPLACED} ${ARGN})
endmacro(blocktree_build_test)

include(CTest)

find_package(GTest REQUIRED)
pasta_block_tree_build_test(block_tree/block_tree_fp_test)
pasta_block_tree_build_test(block_tree/block_tree_lpf_test)

################################################################################
