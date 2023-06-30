/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2022 Daniel Meyer
 *
 * pasta::block_tree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pasta::block_tree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::block_tree.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <omp.h>
#include <vector>
#include <stack>
#include <cmath>

const int BLOCK_SIZE = 8192;

inline int32_t get_left(int32_t i) {
    return i << 1;
}

inline int32_t get_right(int32_t i) {
    return (i << 1) | 1;
}

inline int32_t get_parent(int32_t i) {
    return i >> 1;
}

inline int64_t get_left(int64_t i) {
    return i << 1;
}

inline int64_t get_right(int64_t i) {
    return (i << 1) | 1;
}

inline int64_t get_parent(int64_t i) {
    return i >> 1;
}

inline int32_t getLeft_opt(std::vector<std::vector<int32_t>> &table, int32_t depth, int32_t index, int32_t start) {
    int32_t value = table[0][index];
    if (value == table[depth - 1][0]) return -1;

    int32_t cur = get_parent(start), d, dist = 2;
    for (d = 1; d < depth; d++) {
        if ((cur + 1) * dist > index + 1) cur--;
        if (cur < 0) return -1;

        if (table[d][cur] >= value) cur = get_parent(cur);
        else break;

        dist <<= 1;
    }

    for (; d > 0; d--) {
        if (table[d - 1][get_right(cur)] < value) cur = get_right(cur);
        else cur = get_left(cur);
    }
    return cur;
}

inline int64_t getLeft_opt(std::vector<std::vector<int64_t>> &table, int64_t depth, int64_t index, int64_t start) {
    int64_t value = table[0][index];
    if (value == table[depth - 1][0]) return -1;

    int64_t cur = get_parent(start), d, dist = 2;
    for (d = 1; d < depth; d++) {
        if ((cur + 1) * dist > index + 1) cur--;
        if (cur < 0) return -1;

        if (table[d][cur] >= value) cur = get_parent(cur);
        else break;

        dist <<= 1;
    }

    for (; d > 0; d--) {
        if (table[d - 1][get_right(cur)] < value) cur = get_right(cur);
        else cur = get_left(cur);
    }
    return cur;
}

inline int64_t getRight_opt(std::vector<std::vector<int64_t>> &table, int64_t depth, int64_t index, int64_t start) {
    int64_t value = table[0][index];
    if (value == table[depth - 1][0]) return -1;

    int64_t cur = get_parent(start), d, dist = 2;
    for (d = 1; d < depth; d++) {
        if (cur * dist < index) cur++;
        if (cur * dist >= table[0].size()) return -1;

        if (table[d][cur] >= value) cur = get_parent(cur);
        else break;

        dist <<= 1;
    }

    for (; d > 0; d--) {
        if (table[d - 1][get_left(cur)] < value) cur = get_left(cur);
        else cur = get_right(cur);
    }
    return cur;
}

inline int32_t getRight_opt(std::vector<std::vector<int32_t>> &table, int32_t depth, int32_t index, int32_t start) {
    int32_t value = table[0][index];
    if (value == table[depth - 1][0]) return -1;

    int32_t cur = get_parent(start), d, dist = 2;
    for (d = 1; d < depth; d++) {
        if (cur * dist < index) cur++;
        if (cur * dist >= table[0].size()) return -1;

        if (table[d][cur] >= value) cur = get_parent(cur);
        else break;

        dist <<= 1;
    }

    for (; d > 0; d--) {
        if (table[d - 1][get_left(cur)] < value) cur = get_left(cur);
        else cur = get_right(cur);
    }
    return cur;
}

template<typename size_type>
int32_t
ansv(std::vector<size_type> &array, std::vector<size_type> &ln, std::vector<size_type> &rn, size_type offset,
     size_type n) {
    std::stack<size_type> stack_left;
    for (size_type i = 0; i < n; i++) {
        while (!stack_left.empty() && array[stack_left.top() + offset] > array[i + offset]) {
            stack_left.pop();
        }
        if (stack_left.empty()) ln[i + offset] = -1;
        else ln[i + offset] = stack_left.top() + offset;
        stack_left.push(i);
    }
    std::stack<size_type> stack_right;
    for (size_type i = n - 1; i >= 0; i--) {
        while (!stack_right.empty() && array[stack_right.top() + offset] > array[i + offset]) {
            stack_right.pop();
        }
        if (stack_right.empty()) {
            rn[i + offset] = -1;
        } else {
            rn[i + offset] = stack_right.top() + offset;
        }
        stack_right.push(i);
    }
    return 0;
}

template<typename size_type>
void ansv_omp(std::vector<size_type> &array, std::vector<size_type> &left, std::vector<size_type> &right,
              size_type threads) {
    size_type l2 = std::ceil(log2(array.size()));
    size_type depth = l2 + 1;
    std::vector<std::vector<size_type>> table(depth);
    table[0] = array;
    auto m = (size_type) array.size();
    for (size_type i = 1; i < depth; i++) {
        m = (m + 1) / 2;
        table[i] = std::vector<size_type>(m);
    }
    m = (size_type) array.size();
    for (size_type d = 1; d < depth; d++) {
        size_type m2 = m / 2;
        omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(m2, table, d)
        for (size_type i = 0; i < m2; i++) {
            table[d][i] = std::min(table[d - 1][get_left(i)], table[d - 1][get_right(i)]);
        }

        if (m % 2) {
            table[d][m2] = table[d - 1][get_left(m2)];
        }

        m = (m + 1) / 2;
    }
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(table, array, left, right, depth)
    for (size_type i = 0; i < array.size(); i += BLOCK_SIZE) {
        size_type j = std::min(i + BLOCK_SIZE, (size_type) array.size());
        ansv(array, left, right, i, j - i);

        size_type tmp = i;
        for (size_type k = i; k < j; k++) {
            if (left[k] == -1) {
                if (tmp != -1 && array[tmp] >= array[k]) {
                    tmp = getLeft_opt(table, depth, k, tmp);
                }
                left[k] = tmp;
            }
        }

        tmp = j - 1;
        for (size_type k = j - 1; k >= i; k--) {
            if (right[k] == -1) {
                if (tmp != -1 && array[tmp] >= array[k]) {
                    tmp = getRight_opt(table, depth, k, tmp);
                }
                right[k] = tmp;
            }
        }
    }
}

/******************************************************************************/
