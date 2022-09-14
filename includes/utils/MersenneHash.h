#include <cstddef>
#include <iostream>
#include <vector>
#ifndef MERSENNE_HASH_H
#define MERSENNE_HASH_H

template <typename T>
class MersenneHash
{
    public:
        std::vector<T> const& text_;
        uint64_t hash_;
        uint32_t start_;
        uint32_t length_;
        MersenneHash(std::vector<T> const& text, size_t hash, uint64_t start, uint64_t length) : text_(text), hash_(hash), start_(start), length_(length) {};
        bool operator==(const MersenneHash &other ) const {
//            std::cout << start_ << " " << other.start_ << std::endl;
            if (length_ != other.length_) return false;
            
            for (int i = 0; i < length_; i++) {
                if (text_[start_ + i] != other.text_[other.start_ + i]) {
                    return false;
                }
            }
            return true;
        }
};


namespace std {
    template<typename T>
    struct hash<MersenneHash<T>> {
        std::size_t operator()(const MersenneHash<T>& hS) const {
            return hS.hash_;
        }
    };
}
#endif // MERSENNE_HASH_H