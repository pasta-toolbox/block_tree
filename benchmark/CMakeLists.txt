set(benches
        foo
        bench_lz
        bench_fp
        )
foreach(bench ${benches})
    add_executable(${bench} ${bench}.cpp)
    target_link_libraries(${bench} pasta_bit_vector)
    target_link_libraries(${bench} sdsl)
    target_link_libraries(${bench} tlx)
    target_link_libraries(${bench} sais)
    target_link_libraries(${bench} sais64)
    target_link_libraries(${bench} ${CMAKE_SOURCE_DIR}/extlib/malloc_count/malloc_count.o)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/bit_vector/include/pasta/bit_vector)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/tlx/tlx)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/bit_vector/include/pasta/bit_vector/support)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/sdsl/include/sdsl)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/includes/utils)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/includes/bv_block_tree)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/libsais_cmake/src)
    target_include_directories(${bench} PUBLIC ${CMAKE_SOURCE_DIR}/extlib/malloc_count)
    endforeach()