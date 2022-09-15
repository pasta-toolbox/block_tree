
#include <libsais.h>
#include <libsais64.h>
#include <iostream>
#include <vector>
#include <chrono>
#ifndef BLOCK_TREE_LPFARRAY_H
#define BLOCK_TREE_LPFARRAY_H


template<typename size_type>
int32_t calculate_lz_factor(size_type &z, std::vector<size_type> &lpf, std::vector<size_type> &lz) {
    size_type i = 0;
    lz.push_back(0);
    while (lz[i] < lpf.size() - 1) {
        lz.push_back(lz[i] + std::max(1, lpf[lz[i]]));
        i++;
    }
    z = i;
//    std::cout << "Given Text has " << z << " LZ-factors"<< std::endl;
    return 0;

};

int32_t lpf_array(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &lpf_ptr) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> plcp(text.size());
    std::vector<int64_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais64(text.data(), sa.data() , text.size(), 0, NULL);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais64_plcp(text.data(), sa.data(), plcp.data(), text.size());
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais64_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    auto t04 = std::chrono::high_resolution_clock::now();
    auto ms_int4 = std::chrono::duration_cast<std::chrono::milliseconds>(t04 - t03);
//    std::cout << "LCP TIME " << ms_int4.count() << std::endl;
    std::vector<int64_t> isu(text.size());
    std::vector<int64_t> prev(text.size());
    std::vector<int64_t> next(text.size());
    for(int64_t i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int64_t r = 0; r < text.size() - 1; r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    auto t05 = std::chrono::high_resolution_clock::now();
    auto ms_int5 = std::chrono::duration_cast<std::chrono::milliseconds>(t05 - t04);
//    std::cout << "ALLOCATION TIME " << ms_int5.count() << std::endl;
    for (int64_t i = text.size() - 1; i >= 0; i--) {
        int64_t r = isu[i];
        int64_t next_r = next[r];
        int64_t prev_r = prev[r];
        int64_t lcp_r = lcp[r];
        int64_t lcp_nr = lcp[next_r];
        if (lcp_r <= lcp_nr) {
            lpf[i] = lcp_nr;
            lpf_ptr[i] = sa[next_r];
            lcp[next_r] = lcp_r;
        } else {
            lpf[i] = lcp_r;
            lcp[next_r] = lcp_nr;
            lpf_ptr[i] = sa[prev_r];
        }
        if (prev_r >= 0) {
            next[prev_r] = next_r;
        }
        if (next_r < text.size()) {
            prev[next_r] = prev_r;
        }
    }
    auto t06 = std::chrono::high_resolution_clock::now();
    auto ms_int6 = std::chrono::duration_cast<std::chrono::milliseconds>(t06 - t05);
//    std::cout << "LPF TIME " << ms_int6.count() << std::endl;
    return 0;
}

int32_t lpf_array(std::vector<uint8_t> &text, std::vector<int32_t> &lpf, std::vector<int32_t> &lpf_ptr) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais(text.data(), sa.data() , text.size(), 0, NULL);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais_plcp(text.data(), sa.data(), plcp.data(), text.size());
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    auto t04 = std::chrono::high_resolution_clock::now();
    auto ms_int4 = std::chrono::duration_cast<std::chrono::milliseconds>(t04 - t03);
//    std::cout << "LCP TIME " << ms_int4.count() << std::endl;
    std::vector<int32_t> isu(text.size());
    std::vector<int32_t> prev(text.size());
    std::vector<int32_t> next(text.size());
    for(int32_t i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int32_t r = 0; r < text.size() - 1; r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    auto t05 = std::chrono::high_resolution_clock::now();
    auto ms_int5 = std::chrono::duration_cast<std::chrono::milliseconds>(t05 - t04);
//    std::cout << "ALLOCATION TIME " << ms_int5.count() << std::endl;
    for (int32_t i = text.size() - 1; i >= 0; i--) {
        int32_t r = isu[i];
        int32_t next_r = next[r];
        int32_t prev_r = prev[r];
        int32_t lcp_r = lcp[r];
        int32_t lcp_nr = lcp[next_r];
        if (lcp_r <= lcp_nr) {
            lpf[i] = lcp_nr;
            lpf_ptr[i] = sa[next_r];
            lcp[next_r] = lcp_r;
        } else {
            lpf[i] = lcp_r;
            lcp[next_r] = lcp_nr;
            lpf_ptr[i] = sa[prev_r];
        }
        if (prev_r >= 0) {
            next[prev[r]] = next[r];
        }
        if (next_r < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    auto t06 = std::chrono::high_resolution_clock::now();
    auto ms_int6 = std::chrono::duration_cast<std::chrono::milliseconds>(t06 - t05);
//    std::cout << "LPF TIME " << ms_int6.count() << std::endl;
    return 0;
}
int32_t lpf_array(std::string &text, std::vector<int32_t> &lpf, std::vector<int32_t> &lpf_ptr) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    libsais(reinterpret_cast<const uint8_t *>(text.c_str()), sa.data() , text.size(), 0, NULL);
    libsais_plcp(reinterpret_cast<const uint8_t *>(text.c_str()), sa.data(), plcp.data(), text.size());
    libsais_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    std::vector<int32_t> isu(text.size());
    std::vector<int32_t> prev(text.size());
    std::vector<int32_t> next(text.size());
    for(int i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int r = 0; r < text.size() - 1; r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    for (int i = text.size() - 1; i >= 0; i--) {
        int r = isu[i];
        lpf[i] = std::max(lcp[r], lcp[next[r]]);
        lcp[next[r]] = std::min(lcp[r], lcp[next[r]]);
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    return 0;
}
#endif //BLOCK_TREE_LPFARRAY_H
