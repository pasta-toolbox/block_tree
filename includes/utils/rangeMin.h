//
// Created by daniel on 13.11.22.
//
#include <vector>
#include <math.h>
#include <omp.h>
#ifndef BLOCK_TREE_RANGEMIN_H
#define BLOCK_TREE_RANGEMIN_H
template<typename size_type>
class min_range_q {
public:
    std::vector<size_type>& array;
    size_type m;
    size_type n;
    size_type** table;
    size_type depth;
    size_type BSIZE = 16;
    size_type threads = 1;
    explicit min_range_q(std::vector<size_type>& _array, size_type _n, size_type _threads) : array(_array), n(_n), threads(_threads){
        m = 1 + (array.size() - 1)/BSIZE;
        precompute();
    };
    void precompute() {
        depth = log2(m) + 1;
        table = new size_type*[depth];
        omp_set_num_threads(threads);
#pragma omp parallel for default(none)
            for (size_type i = 0; i < depth; i++) {
                table[i] = new size_type[n];
            }


#pragma omp parallel for default(none)
            for (size_type i = 0; i < m; i++) {
                size_type start = i * BSIZE;
                size_type end = std::min(start + BSIZE, n);
                size_type k = i * BSIZE;
                for (size_type j = start + 1; j < end; j++)
                    if (array[j] < array[k]) k = j;
                table[0][i] = k;
            }

        size_type dist = 1;
        for(size_type j=1;j<depth;j++) {

#pragma omp parallel for default(none) shared(dist, j)
            for (size_type i = 0; i < m - dist; i++) {
                if (array[table[j - 1][i]] <= array[table[j - 1][i + dist]])
                    table[j][i] = table[j - 1][i];
                else table[j][i] = table[j - 1][i + dist];
            }


#pragma omp parallel for default(none) shared(dist, j)
            for (size_type i = m - dist; i < m; i++) {
                table[j][i] = table[j - 1][i];
            }

            dist *= 2;
        }
    };
    size_type query(size_type i, size_type j) {
            //same block
            if (j-i < BSIZE) {
                size_type r = i;
                for (size_type k = i+1; k <= j; k++)
                    if (array[k] < array[r]) r = k;
                return r;
            }
        size_type block_i = i/BSIZE;
        size_type block_j = j/BSIZE;
        size_type min = i;
            for(size_type k=i+1;k<(block_i+1)*BSIZE;k++){
                if(array[k] < array[min]) min = k;
            }
            for(size_type k=j; k>=(block_j)*BSIZE;k--){
                if(array[k] < array[min]) min = k;
            }
            if(block_j == block_i + 1) return min;
        size_type outOfBlockMin;
            //not same or adjacent blocks
            if(block_j > block_i + 1){
                block_i++;
                block_j--;
                if(block_j == block_i) outOfBlockMin = table[0][block_i];
                else if(block_j == block_i + 1) outOfBlockMin = table[1][block_i];
                else {
                    size_type k = log2(block_j - block_i);
                    size_type p = 1<<k; //2^k
                    outOfBlockMin = array[table[k][block_i]] <= array[table[k][block_j+1-p]]
                                    ? table[k][block_i] : table[k][block_j+1-p];
                }
            }

            return array[min] < array[outOfBlockMin] ? min : outOfBlockMin;
    }
    ~min_range_q() {
        for(size_type i=0;i<depth;i++){
            delete[] table[i];
        }
        delete[] table;
    };
};



#endif //BLOCK_TREE_RANGEMIN_H
