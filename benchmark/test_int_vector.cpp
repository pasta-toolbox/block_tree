//
// Created by daniel on 20.09.22.
//
#include <sdsl/int_vector.hpp>
#include <vector>
#include <unordered_map>
int main() {
    auto vector = std::vector<std::vector<sdsl::int_vector<>*>*>();
    vector.resize(6, new std::vector<sdsl::int_vector<>*>());
    std::cout << vector.size() << std::endl;
    for (auto a: vector) {
        a->resize(4);
        std::cout << a->size() << std::endl;
        for (auto aa: *a) {
            aa = new sdsl::int_vector<>(5,0);
            std::cout << aa->size() << std::endl;
        }

    }
    return 0;
}
