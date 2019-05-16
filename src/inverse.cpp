//
// Created by Paolo on 16/05/2019.
//

#include <cassert>
#include "inverse.h"
#include "resol.h"

Vector e(size_t i, size_t n) {
    assert(i < n);
    Vector e(0.0, n);
    e[i] = 1.0;
    return e;
}

Matrix identity(size_t n) {
    Matrix Id(n);
    for (index_t i = 0; i < n; ++i) Id(i, i) = 1.0;
    return Id;
}

Matrix inverse(const LUDecomposition &luObj) {
    const size_t n = luObj.decompMatrix().size();
    Matrix Inv(n);
    for (index_t i = 0; i < n; ++i) {
        Vector column = solveLU(luObj, e(i, n));
        for (index_t j = 0; j < n; ++j) Inv[j][i] = column[j];
    }
    return Inv;
}

Matrix inverse(const Matrix &A) {
    return inverse(LUDecomposition(A));
}
