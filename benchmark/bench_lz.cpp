//
// Created by daniel on 03.09.22.
//
#include <iostream>
#include "libsais.h"
#include "pasta/block_tree/utils/lpf_array.hpp"
#include <vector>
#include <pasta/bit_vector/bit_vector.hpp>
#include <tlx/cmdline_parser.hpp>
#include <fstream>
#include <sstream>
#include "pasta/block_tree/construction/block_tree_lpf.hpp"
#include <chrono>
#include <type_traits>
#include <iostream>

int main(int argc, char* argv[]) {
    tlx::CmdlineParser cp;
    // add a byte size argument which the user can enter like '1gi'
    uint64_t a_size = 0;
    cp.add_bytes('s', "size", a_size,
                 "Number of bytes to process.");

    // process command line
    if (!cp.process(argc, argv))
        return -1; // some error occurred and help was always written to user.

//    std::cout << "Command line parsed okay." << std::endl;
    std::cout <<  std::endl  << "Run with "<< a_size << " Bytes" << std::endl;
    std::string test(a_size, ' ');
//    std::ifstream t("/home/daniel/blocktree-experiments/data/Escherichia_Coli");
    std::ifstream t("/home/daniel/blocktree-experiments/data/english.1024MB");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.de.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/einstein.en.txt");
//    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
//    std::ifstream t("/Users/daniel/Downloads/einstein.en.txt");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
//   test = "ABCDABCDEFGHABCDDE12";
    std::vector<uint8_t> vec(test.begin(), test.end());
//    std::vector<int32_t> lpf(test.size());
//    std::vector<int32_t> lpf_ptr(test.size());
//    int64_t lz = 0;
//    int64_t k = 0;
//    while (lz < test.size()) {
//        lz = lz + std::max(1LL, lpf[lz + 1]);
//        k++;
//    }
//    std::cout << "lz " << lz.size() << std::endl;
    std::vector<int32_t> lpf(vec.size());
    std::vector<int32_t> lpf_ptr(vec.size());
    std::vector<int32_t> lz;
    int32_t lzn = 0;
    lpf_array(vec, lpf, lpf_ptr);
    calculate_lz_factor(lzn,lpf, lz);
    auto t01 = std::chrono::high_resolution_clock::now();
    BlockTreeLPF<uint8_t, int32_t>*  bt = new BlockTreeLPF<uint8_t, int32_t>(vec, 2, 1,1,lpf, lpf_ptr, lz, true,true);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);



    // BV_BlockTree_lpf_heuristic<uint8_t, int32_t>*  bt2 = new BV_BlockTree_lpf_heuristic<uint8_t, int32_t>(vec, 2, 1,lpf, lpf_ptr, lz);
    // auto t0x = std::chrono::high_resolution_clock::now();
    // BV_BlockTree_lpf_theory<uint8_t, int32_t>*  bt3 = new BV_BlockTree_lpf_theory<uint8_t, int32_t>(vec, 2, 1,1,lpf, lpf_ptr, lz);

    // auto t0y = std::chrono::high_resolution_clock::now();
    // auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t0x - t02);
    // auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t0y - t0x);
    // std::cout << "time pruned " <<  ms_int.count() << std::endl;
    // std::cout << "time appr " <<  ms_int2.count() << std::endl;
    // std::cout << "time theory " <<  ms_int3.count() << std::endl;
//    std::cout << "pruned " << bt->block_size_lvl_[0]  << std::endl;
//    for (auto bv: bt->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
    for (auto bv: bt->block_tree_pointers_) {
        for (auto a: *bv ) {
            std::cout << a << " ";
        }
        std::cout << std::endl << std::endl;
    }
//    for (auto bv: bt->block_tree_offsets_) {
//        for (auto a: *bv ) {
//            std::cout << a << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }
//    std::cout << "theo " << bt3->block_size_lvl_[0] <<  std::endl;
//    for (auto bv: bt3->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
//    for (auto bv: bt3->block_tree_pointers_) {
//        for (auto a: *bv ) {
//            std::cout << a << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }
//    for (auto bv: bt3->block_tree_offsets_) {
//        for (auto a: *bv ) {
//            std::cout << a << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }
//    std::cout << "appr "<< bt2->block_size_lvl_[0] << std::endl;
//    for (auto bv: bt2->block_tree_types_) {
//        std::cout << *bv << std::endl;
//    }
//    for (auto bv: bt2->block_tree_pointers_) {
//        for (auto a: *bv ) {
//            std::cout << a << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }
//    for (auto bv: bt2->block_tree_offsets_) {
//        for (auto a: *bv ) {
//            std::cout << a << " ";
//        }
//        std::cout << std::endl << std::endl;
//    }
//    for (auto l: bt->leaves_) {
//        std::cout << l;
//    }
//    std::cout << std::endl;
//    for (auto l: bt2->leaves_) {
//        std::cout << l;
//    }
//    std::cout << std::endl;
    int j = 0;
//    std::cout << 146 << " " << bt2->access(146) << std::endl;
//    std::cout << 146 << " " << bt->access(146) << std::endl;
    for (uint64_t i = 0; i < test.size(); i++) {
        if (bt->access(i) != vec[i]) {
            std::cout << i << " " << bt->access(i) << " " << vec[i] << std::endl;
            j++;
        }
    }
//
//
//    }
//    std::cout << j << " errors in appr." << std::endl;
//    for (int i = 0; i < test.size(); i++) {
//        if (bt3->access(i) != vec[i]) {
//            std::cout << bt3->access(i) << " " <<  vec[i] << std::endl;
//            j++;
//        }
//    }
//    std::cout << j << " errors in theo." << std::endl;
//
//
//    auto sum = 0;
////    for (auto a: bt->block_tree_types_) {
////        std::cout << *a << std::endl;
////    }
////    std::cout << "=======" << std::endl;
////    for (auto a: bt2->block_tree_types_) {
////        std::cout << *a << std::endl;
////    }
//    auto t03 = std::chrono::high_resolution_clock::now();
//
//    for (int i = 0; i < test.size(); i++) {
//        auto x = bt->access(i);
//        if (x != vec[i]) {
//            j++;
//            std::cout << i<< " "<< x << " " <<  vec[i] << std::endl;
//        }
//
//        sum += x;
//    }
//    std::cout << j << " errors in pruned." << std::endl;
//    auto t04 = std::chrono::high_resolution_clock::now();
//    auto ns_int2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t04 - t03);
//
//    std::cout << "Errors: " <<  j << " avg time: " << (double) (ns_int2.count()/test.size()) << std::endl;
//    for (auto v: bt->block_tree_pointers_) {
//        std::cout << (int)v->width() << std::endl;
//    }
    std::cout << "pruned" << std::endl;
    for (uint64_t i = 0; i < bt->block_tree_types_.size(); i++) {
        auto& lvl_bv = *bt->block_tree_types_[i];
        auto& lvl_bv_rs = *bt->block_tree_types_rs_[i];
        auto& lvl_ptr = *bt->block_tree_pointers_[i];
        auto& lvl_off = *bt->block_tree_offsets_[i];
        std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
        std::cout << lvl_bv.size() << "/";
        std::cout << lvl_bv_rs.space_usage() << "/";
        std::cout << "("<<lvl_ptr.size()<<","<<(int) lvl_ptr.width()<< ","<<lvl_ptr.bit_size()<<")/";
        std::cout << "("<<lvl_off.size()<<","<< (int) lvl_off.width()<< ","<<lvl_off.bit_size()<<")"<<std::endl;
    }
    std::cout << " leaves: " << bt->leaves_.size() << std::endl;
    std::cout << "appr" << std::endl;
    // for (int i = 0; i < bt2->block_tree_types_.size(); i++) {
    //     auto& lvl_bv = *bt2->block_tree_types_[i];
    //     auto& lvl_bv_rs = *bt2->block_tree_types_rs_[i];
    //     auto& lvl_ptr = *bt2->block_tree_pointers_[i];
    //     auto& lvl_off = *bt2->block_tree_offsets_[i];
    //     std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
    //     std::cout << lvl_bv.size() << "/";
    //     std::cout << lvl_bv_rs.space_usage() << "/";
    //     std::cout << "("<<lvl_ptr.size()<<","<<(int) lvl_ptr.width()<< ","<<lvl_ptr.bit_size()<<")/";
    //     std::cout << "("<<lvl_off.size()<<","<< (int) lvl_off.width()<< ","<<lvl_off.bit_size()<<")"<<std::endl;
    // }
    // std::cout << " leaves: " << bt2->leaves_.size() << std::endl;
    // std::cout << "theo" << std::endl;
    // for (int i = 0; i < bt3->block_tree_types_.size(); i++) {
    //     auto& lvl_bv = *bt3->block_tree_types_[i];
    //     auto& lvl_bv_rs = *bt3->block_tree_types_rs_[i];
    //     auto& lvl_ptr = *bt3->block_tree_pointers_[i];
    //     auto& lvl_off = *bt3->block_tree_offsets_[i];
    //     std::cout << i << " lvl/bv_s/bv_rs/pointer(n,w,s)/pointer(n,w,s) ";
    //     std::cout << lvl_bv.size() << "/";
    //     std::cout << lvl_bv_rs.space_usage() << "/";
    //     std::cout << "("<<lvl_ptr.size()<<","<<(int) lvl_ptr.width()<< ","<<lvl_ptr.bit_size()<<")/";
    //     std::cout << "("<<lvl_off.size()<<","<< (int) lvl_off.width()<< ","<<lvl_off.bit_size()<<")"<<std::endl;
    // }
    // std::cout << " leaves: " << bt3->leaves_.size() << std::endl;
    bt->print_space_usage();
    // bt2->print_space_usage();
    // bt3->print_space_usage();
    // int sizer = 0;
    // int sizer2 = 0;
    // for (auto pv: bt->block_tree_pointers_) {
    //     auto& ptr = *pv;
    //     auto x = 0;
    //     for (auto p : *pv) {
    //         if (p > x) {
    //             x = p;
    //         }
    //     }
    //     std::cout << (int)ptr.width() << " " << 32 -  __builtin_clz(static_cast<unsigned int>(x) | 1) << std::endl;
    //     sizer +=  (32 -  __builtin_clz(static_cast<unsigned int>(x) | 1)) * ptr.size();
    //     sizer2 += (int) ptr.width() * ptr.size();
    // }
    // std::cout << sizer << " " << sizer2 << " " << sizer2 - sizer << std::endl;
    // auto x = 0;
    delete bt;
    // auto y = 0;
    // delete bt2;
    // auto z = 0;
    // delete bt3;
    // auto a = 0;
    // std::cout <<"Pruned: " <<  (double )8 * (x - y)/test.size() << " Bits/S" << std::endl;
    // std::cout <<"Appr.: " <<  (double )8 * (y - z)/test.size() << " Bits/S" << std::endl;
    // std::cout <<"Theory: " <<  (double )8 * (z - a)/test.size() << " Bits/S" << std::endl;
//    std::cout << test << std::endl;

    return 0;
}
