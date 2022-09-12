#include <iostream>

#ifndef MERSENNE_RABINKARP_H
#define MERSENNE_RABINKARP_H

template <class T, class size_type>
class MersenneRabinKarp
{
    __extension__ typedef unsigned __int128 uint128_t;
public:
    uint64_t hash_;
    static constexpr uint128_t  prime_ = 2305843009213693951ULL;
    uint128_t sigma_;
    uint128_t max_sigma_;
    uint64_t init_;
    uint64_t length_;
    std::vector<T> const& text_;
    MersenneRabinKarp(std::vector<T> const& text, uint64_t sigma, uint64_t init, uint64_t length) : text_(text), sigma_(sigma), init_(init), length_(length)
    {
        max_sigma_ = 1;
        uint128_t fp = 0;
        uint128_t sigma_c = 1;
        static constexpr uint128_t kPrime = prime_;
        for (int i = init_ ; i <  init_ + length_; i++)
        {
            fp = mersenneModulo(fp * sigma);
            fp = mersenneModulo(fp + text_[i]);
        }
        for (int i = 0; i < length_ - 1; i++) {
            sigma_c =  mersenneModulo(sigma_c * sigma_);
        }
        hash_ = static_cast<uint64_t>(fp);
        max_sigma_ = static_cast<uint64_t>(sigma_c);
    };
    void restart(uint64_t index) {
        if (index >= text_.size()) {
            return;
        }
        init_ = index;
        max_sigma_ = 1;
        uint128_t fp = 0;
        uint128_t sigma_c = 1;
        static constexpr uint128_t kPrime = prime_;
        for (int i = init_ ; i <  init_ + length_; i++)
        {
            fp = mersenneModulo(fp * sigma_);
            fp = mersenneModulo(fp + text_[i]);
        }
        for (int i = 0; i < length_ - 1; i++) {
            sigma_c =  mersenneModulo(sigma_c * sigma_);
        }
        hash_ = static_cast<uint64_t>(fp);
        max_sigma_ = static_cast<uint64_t>(sigma_c);
    };
    inline uint128_t mersenneModulo(uint128_t k) {
        return k % prime_;
//        uint128_t i = (k & prime_) + (k >> power_);
//        return (i >= prime_) ? i - prime_ : i;
    };
    void next() {
        if (text_.size() == init_+1) {
            return;
        }
        uint128_t fp = hash_;
        T out_char = text_[init_];
        T in_char = text_[init_ + length_];
        uint128_t out_char_influence = out_char * max_sigma_;
        out_char_influence = mersenneModulo(out_char_influence);
        if(out_char_influence < hash_) {
            fp -= out_char_influence;
        } else {
            fp = prime_ - (out_char_influence - fp);
        }
        fp *= sigma_;
        fp += in_char;
        fp = mersenneModulo(fp);
        hash_ = static_cast<uint64_t>(fp);
        init_ ++;
    };
    uint64_t hash() {
        return hash_;
    };
    std::string toString(__uint128_t num) {
        std::string str;
        do {
            int digit = num % 10;
            str = std::to_string(digit) + str;
            num = (num - digit) / 10;
        } while (num != 0);
        return str;
    }
};

#endif // MERSENNE_RABINKARP_H