#include <iostream>
#include <libsais.h>
#include <lpfarray.hpp>
#include <vector>
#include <bit_vector.hpp>
#include <tlx/cmdline_parser.hpp>
#include <fstream>
#include <sstream>
#include <bv_blocktree_lpf.hpp>
#include "malloc_count/malloc_count.h"
#include <chrono>
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
    std::ifstream t("/home/daniel/blocktree-experiments/data/influenza");
    std::stringstream buffer;
    t.read(&test[0], a_size);
//    test = "NNBOBOTWNNBOBIOOTBSHTFNEBOBOTWNEBOBOTWNEBOBIOOTBSHTFNSBOBOTW";
    std::vector<uint8_t> vec(test.begin(), test.end());
    std::vector<int64_t> lpf(test.size());
    std::vector<int64_t> lpf_ptr(test.size());
//    int64_t lz = 0;
//    int64_t k = 0;
//    while (lz < test.size()) {
//        lz = lz + std::max(1LL, lpf[lz + 1]);
//        k++;
//    }
//    std::cout << "lz " << lz.size() << std::endl;
    auto t01 = std::chrono::high_resolution_clock::now();
    BVBlockTree bt(vec, 2, 32,1);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "time " <<  ms_int.count() << std::endl;

//    int childs = 1;
//    int blocks = 1;
//    int child_size = test.size()/childs;
//    std::vector<int64_t> block_text_inx(childs);
//    for (int i = 0; i < childs; i++) {
//        block_text_inx[i] = (i * child_size);
//    }
//    std::vector<pasta::BitVector*> block_types;
//    std::vector<std::vector<int64_t>*> block_pointers;
//    std::vector<std::vector<int64_t>*> block_offsets;
//    while (child_size > 1) {
//        std::cout << child_size << std::endl;
//        pasta::BitVector* bv = new pasta::BitVector(block_text_inx.size(),0);
//        std::vector<int64_t>* pointers = new std::vector<int64_t>();
//        std::vector<int64_t>* offsets = new std::vector<int64_t>();
//        std::vector<int64_t> block_text_inx_new;
//        for (auto a : block_text_inx) {
//            std::cout << a <<",";
//        }
//        std::cout << std::endl;
//        for(int i = block_text_inx.size() - 1; i >=  0; i--) {
//            // if block is already marked we don't have to check for replacements
//            if ((*bv)[i] == 1) {
//                continue;
//            }
//            // ind is initially block i's starting pos in the text
//            int ind = block_text_inx[i];
//            bool has_ptr = false;
//            // we look until we find leftmost occ of [ind..ind+blocksize]
//            while (lpf[ind] >= child_size) {
//                if (child_size + lpf_ptr[ind] >= ind) {
//                    // our found occs overlapps with the block
//                    // we will try if we find another one
//                    ind = lpf_ptr[ind];
//                } else if (lpf[lpf_ptr[ind]] < child_size) {
//                    // we found the leftmost occ
//                    has_ptr = true;
//                    // determine which block(s) it points at
//                    int ptr = lpf_ptr[ind];
//                    for (int j = 1; j < block_text_inx.size(); j++) {
//                        if (ptr < block_text_inx[j]) {
//                            //std::cout << "found a block " << i << " pointing to Block " << j-1 <<  std::endl;
//                            (*bv)[j-1] = 1;
//                            if (ptr % child_size != 0) {
//                                (*bv)[j] = 1;
//                            }
//                            int64_t p = j - 1;
//                            pointers->push_back(p);
//                            offsets->push_back(ptr % child_size);
//                            break;
//                        }
//                    }
//                    break;
//                } else {
//                    ind = lpf_ptr[ind];
//                }
//            }
//            // if there is no block its pointing at it is its leftmost occ
//            if (!has_ptr) {
//                (*bv)[i] = 1;
//            }
//        }
//        for (int i = 0; i < bv->size(); i++) {
//            if ((*bv)[i] == 1) {
//                for (int j = 0; j < 2; j++) {
//                    block_text_inx_new.push_back(block_text_inx[i] + j * child_size/2);
//                }
//            }
//        }
//        std::cout << "HI " << bv->size()  << std::endl;
//        block_types.push_back(bv);
//        block_pointers.push_back(pointers);
//        block_offsets.push_back(offsets);
//        block_text_inx = block_text_inx_new;
//        child_size = child_size/2;
//        std::cout << child_size << std::endl;
//    }
//    std::cout << "lvl _sizes" << std::endl;
//    for (auto a: block_types) {
//        std::cout << a->size() << std::endl;
//    }
//    std::cout << "offis" << std::endl;
//    for (auto a : block_offsets) {
//        for (int i = 0; i < a->size(); i++) {
//            std::cout << (*a)[i] << ",";
//        }
//        std::cout << std::endl;
//    }
    return 0;
}