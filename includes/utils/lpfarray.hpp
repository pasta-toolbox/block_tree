
#include <libsais.h>
#include <libsais64.h>
#include <stack>
#include <iostream>
#include <vector>
#include <chrono>
#include <rangeMin.h>
#include <cmath>
#include <omp.h>
#include <ANSV.h>
#ifndef BLOCK_TREE_LPFARRAY_H
#define BLOCK_TREE_LPFARRAY_H


template<typename size_type>
int32_t calculate_lz_factor(size_type &z, std::vector<size_type> &lpf, std::vector<size_type> &lz) {
    size_type i = 0;
    lz.push_back(0);
    while (lz[i] < lpf.size() - 1) {
        lz.push_back(lz[i] + std::max(1L, (int64_t) lpf[lz[i]]));
        i++;
    }
    z = i;
//    std::cout << "Given Text has " << z << " LZ-factors"<< std::endl;
    return 0;

};
template<typename size_type>
int32_t calculate_number_lz_factor(size_type &z, std::vector<size_type> &lpf) {
    size_type counter = 0;
    size_type lz = 0;
    while (lz < lpf.size() - 1) {
        lz += std::max(1, lpf[lz]);
        counter++;
    }
    z = counter;
//    std::cout << "Given Text has " << z << " LZ-factors"<< std::endl;
    return 0;

};
int32_t lpf_array(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &lpf_ptr) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> p_lcp(text.size());
    std::vector<int64_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais64(text.data(), sa.data() , text.size(), 0, NULL);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais64_plcp(text.data(), sa.data(), p_lcp.data(), text.size());
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais64_lcp(p_lcp.data(), sa.data(), lcp.data(), text.size());
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
int32_t lpf_array_omp(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &lpf_ptr, int32_t threads) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> plcp(text.size());
    std::vector<int64_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais64_omp(text.data(), sa.data() , text.size(), 0, NULL, threads);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais64_plcp_omp(text.data(), sa.data(), plcp.data(), text.size(),threads);
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais64_lcp_omp(plcp.data(), sa.data(), lcp.data(), text.size(), threads);
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
int32_t lpf_array_omp(std::vector<uint8_t> &text, std::vector<int32_t> &lpf, std::vector<int32_t> &lpf_ptr, int32_t threads) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais_omp(text.data(), sa.data() , text.size(), 0, NULL, threads);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais_plcp_omp(text.data(), sa.data(), plcp.data(), text.size(), threads);
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais_lcp_omp(plcp.data(), sa.data(), lcp.data(), text.size(), threads);
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
        prev[r] = r - 1;
        next[r] = r + 1;
    }
    lcp[lcp.size() - 1] = 0;
    auto t05 = std::chrono::high_resolution_clock::now();
    auto ms_int5 = std::chrono::duration_cast<std::chrono::milliseconds>(t05 - t04);
//    std::cout << "ALLOCATION TIME " << ms_int5.count() << std::endl;
    for (int32_t i = text.size(); i >= 0; i--) {
        int32_t r = isu[i];
        if (lcp[r] <= lcp[next[r]]) {
            lpf[i] = lcp[next[r]];
            lpf_ptr[i] = sa[next[r]];
            lcp[next[r]] = lcp[r];
        } else {
            lpf[i] = lcp[r];
            lpf_ptr[i] = sa[prev[r]];
            lcp[next[r]] = lcp[next[r]];

        }
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    auto t06 = std::chrono::high_resolution_clock::now();
    auto ms_int6 = std::chrono::duration_cast<std::chrono::milliseconds>(t06 - t05);
//    std::cout << "LPF TIME " << ms_int6.count() << std::endl;
    return 0;
}
int32_t lpf_array_stack(std::vector<uint8_t> &text, std::vector<int32_t> &lpf, std::vector<int32_t> &lpf_ptr) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais(text.data(), sa.data() , (int) text.size(), 0, nullptr);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais_plcp(text.data(), sa.data(), plcp.data(), (int) text.size());
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais_lcp(plcp.data(), sa.data(), lcp.data(), (int) text.size());
    auto t04 = std::chrono::high_resolution_clock::now();
    std::stack<std::pair<int32_t, int32_t>> stacker;
    sa.push_back(-1);
    lcp.push_back(0);
    stacker.push(std::pair<int32_t, int32_t>(0,sa[0]));
    for (int32_t i = 1; i < sa.size(); i++) {
        int32_t lcp_i = lcp[i];
        while (!stacker.empty() && sa[i] < stacker.top().second) {
            std::pair<int32_t, int32_t> v = stacker.top();
            lpf[v.second] = std::max(v.first, lcp_i);
            lcp_i = std::min(v.first, lcp_i);
            stacker.pop();
            if (lpf[v.second] == 0) {
                lpf_ptr[v.second] = -1;
            } else if (v.first > lcp_i) {
                lpf_ptr[v.second] = stacker.top().second;
            } else {
                lpf_ptr[v.second] = sa[i];
            }

        }
        if (i < sa.size()) {
            stacker.push(std::pair<int32_t, int32_t>(lcp_i, sa[i]));
        }
    }
    return 0;
}
int32_t lpf_array_stack(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &lpf_ptr) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> plcp(text.size());
    std::vector<int64_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais64(text.data(), sa.data() , (int) text.size(), 0, nullptr);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais64_plcp(text.data(), sa.data(), plcp.data(), (int) text.size());
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais64_lcp(plcp.data(), sa.data(), lcp.data(), (int) text.size());
    auto t04 = std::chrono::high_resolution_clock::now();
    std::stack<std::pair<int64_t, int64_t>> stacker;
    sa.push_back(-1);
    lcp.push_back(0);
    stacker.push(std::pair<int64_t, int64_t>(0,sa[0]));
    for (int64_t i = 1; i < sa.size(); i++) {
        int64_t lcp_i = lcp[i];
        while (!stacker.empty() && sa[i] < stacker.top().second) {
            std::pair<int64_t, int64_t> v = stacker.top();
            lpf[v.second] = std::max(v.first, lcp_i);
            lcp_i = std::min(v.first, lcp_i);
            stacker.pop();
            if (lpf[v.second] == 0) {
                lpf_ptr[v.second] = -1;
            } else if (v.first > lcp_i) {
                lpf_ptr[v.second] = stacker.top().second;
            } else {
                lpf_ptr[v.second] = sa[i];
            }

        }
        if (i < sa.size()) {
            stacker.push(std::pair<int64_t, int64_t>(lcp_i, sa[i]));
        }
    }
    lpf_ptr[0] = -1;
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
        prev[r] = r - 1;
        next[r] = r + 1;
    }
    lcp[lcp.size() - 1] = 0;
    auto t05 = std::chrono::high_resolution_clock::now();
    auto ms_int5 = std::chrono::duration_cast<std::chrono::milliseconds>(t05 - t04);
//    std::cout << "ALLOCATION TIME " << ms_int5.count() << std::endl;
    for (int32_t i = text.size(); i >= 0; i--) {
        int32_t r = isu[i];
        if (lcp[r] <= lcp[next[r]]) {
            lpf[i] = lcp[next[r]];
            lpf_ptr[i] = sa[next[r]];
            lcp[next[r]] = lcp[r];
        } else {
            lpf[i] = lcp[r];
            lpf_ptr[i] = sa[prev[r]];
            lcp[next[r]] = lcp[next[r]];

        }
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
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

//int32_t ansv(std::vector<int32_t>& array, std::vector<int32_t>& ln, std::vector<int32_t>& rn, int32_t offset) {
//    std::stack<int32_t> stack_left;
//    for (int32_t i = 0; i < array.size(); i++) {
//        while (!stack_left.empty() && array[stack_left.top()] > array[i]) stack_left.pop();
//        if (stack_left.empty()) ln[i] = -1;
//        else ln[i] = stack_left.top() + offset;
//        stack_left.push(i);
//    }
//    std::stack<int32_t> stack_right;
//    for (int32_t i = (int) array.size() - 1; i >= 0; i--) {
//        while (!stack_right.empty() && array[stack_right.top()] > array[i]) stack_right.pop();
//        if (stack_right.empty()) rn[i] = -1;
//        else rn[i] = stack_right.top() + offset;
//        stack_right.push(i);
//    }
//    return 0;
//}

int32_t lpf_array_ansv(std::vector<uint8_t> &text, std::vector<int32_t> &lpf, std::vector<int32_t> &prev_Occ, int32_t threads) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais_omp(text.data(), sa.data() , text.size(), 0, NULL, threads);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
//    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais_plcp_omp(text.data(), sa.data(), plcp.data(), text.size(), threads);
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
//    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais_lcp_omp(plcp.data(), sa.data(), lcp.data(), text.size(), threads);
    auto t04 = std::chrono::high_resolution_clock::now();
    auto ms_int4 = std::chrono::duration_cast<std::chrono::milliseconds>(t04 - t03);
//    std::cout << "LCP TIME " << ms_int4.count() << std::endl;

    std::vector<int32_t> l(text.size());
    std::vector<int32_t> r(text.size());

    ansv_omp(sa, l, r, threads);

    auto rmq = min_range_q<int32_t>(lcp, lcp.size(), threads);
    omp_set_num_threads(threads);

#pragma omp parallel for default(none) shared(lpf,l,r,sa,lcp,rmq,prev_Occ)
        for (int32_t i = 0; i < lpf.size(); i++) {
            int32_t l_lcp = 0, r_lcp = 0;
            int32_t ln = l[i], rn = r[i];
            int32_t sai = sa[i];
            if (ln != -1) {
                l_lcp = lcp[rmq.query(ln + 1, i)];
            }
            if (rn != -1) {
                r_lcp = lcp[rmq.query(i + 1, rn)];
            }

            if (l_lcp == 0 && r_lcp == 0) {
                prev_Occ[sai] = -1;
                lpf[sai] = 1;
            } else if (l_lcp > r_lcp) {
                prev_Occ[sai] = sa[ln];
                lpf[sai] = l_lcp;
            } else {
                prev_Occ[sai] = sa[rn];
                lpf[sai] = r_lcp;
            }
        }
    return 0;
}
int32_t lpf_array_ansv(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &prev_Occ, int64_t threads) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> plcp(text.size());
    std::vector<int64_t> lcp(text.size());
    auto t01 = std::chrono::high_resolution_clock::now();
    libsais64_omp(text.data(), sa.data() , (int64_t) text.size(), 0, NULL, threads);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto ms_int2 = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << "SA TIME " << ms_int2.count() << std::endl;
    libsais64_plcp_omp(text.data(), sa.data(), plcp.data(), (int64_t) text.size(), threads);
    auto t03 = std::chrono::high_resolution_clock::now();
    auto ms_int3 = std::chrono::duration_cast<std::chrono::milliseconds>(t03 - t02);
    std::cout << "PLCP TIME " << ms_int3.count() << std::endl;
    libsais64_lcp_omp(plcp.data(), sa.data(), lcp.data(), (int64_t) text.size(), threads);
    auto t04 = std::chrono::high_resolution_clock::now();
    auto ms_int4 = std::chrono::duration_cast<std::chrono::milliseconds>(t04 - t03);
    std::cout << "LCP TIME " << ms_int4.count() << std::endl;

    std::vector<int64_t> l(text.size());
    std::vector<int64_t> r(text.size());
    std::cout << "test" << std::endl;
    ansv_omp(sa, l, r, threads);

    auto rmq = min_range_q<int64_t>(lcp, lcp.size(), threads);
    omp_set_num_threads(threads);

#pragma omp parallel for default(none) shared(lpf,l,r,sa,lcp,rmq,prev_Occ)
    for (int32_t i = 0; i < lpf.size(); i++) {
        int32_t l_lcp = 0, r_lcp = 0;
        int32_t ln = l[i], rn = r[i];
        int32_t sai = sa[i];
        if (ln != -1) {
            l_lcp = lcp[rmq.query(ln + 1, i)];
        }
        if (rn != -1) {
            r_lcp = lcp[rmq.query(i + 1, rn)];
        }

        if (l_lcp == 0 && r_lcp == 0) {
            prev_Occ[sai] = -1;
            lpf[sai] = 1;
        } else if (l_lcp > r_lcp) {
            prev_Occ[sai] = sa[ln];
            lpf[sai] = l_lcp;
        } else {
            prev_Occ[sai] = sa[rn];
            lpf[sai] = r_lcp;
        }
    }
    return 0;
}
#endif //BLOCK_TREE_LPFARRAY_H
