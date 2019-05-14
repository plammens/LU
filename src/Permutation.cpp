//
// Created by Paolo on 14/05/2019.
//

#include "Permutation.h"


Permutation::Permutation(size_t n) : _vec(n) {
    for (index_t i = 0; i < n; ++i) _vec[i] = i;
}

void Permutation::permute(index_t a, index_t b) {
    if (a != b) {
        std::swap(_vec[a], _vec[b]);
        _parity = not _parity;
    }
}