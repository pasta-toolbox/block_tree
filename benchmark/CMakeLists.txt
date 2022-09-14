set(benches
        foo
        bench_lz
        bench_fp
        )
foreach(bench ${benches})
    add_executable(${bench} ${bench}.cpp)
    target_link_libraries(${bench} block_tree)
    target_link_libraries(${bench} ${CMAKE_SOURCE_DIR}/extlib/malloc_count/malloc_count.o)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/malloc_count)
    endforeach()