macro(blocktree_build_test TESTNAME)
    string(REPLACE "/" "_" TESTNAME_REPLACED "${TESTNAME}")
    add_executable(${TESTNAME_REPLACED}
            ${TESTNAME}.cpp)
    target_link_libraries(${TESTNAME_REPLACED} PRIVATE
            block_tree
            malloc_count
            GTest::Main)
    include_directories(${TESTNAME_REPLACED} PRIVATE ${gtest_SOURCE_DIR}/include)
    add_test(
            NAME Tester
            COMMAND ${TESTNAME_REPLACED} ${ARGN})
endmacro(blocktree_build_test)

find_package(GTest REQUIRED)
blocktree_build_test(block_tree/block_tree_test)
