//
// Created by Paolo on 17/05/2019.
//

#include <errors.h>
#include "extra.h"

double determinant(const Matrix &A) {
    try { return determinant(LUDecomposition(A)); } 
    catch (SingularMatrixError &) { return 0.0; }
}

double determinant(const LUDecomposition &luObj) {
    const auto &U = luObj.decompMatrix();
    const size_t n = U.size();
    double prod = 1.0;
    for (index_t i = 0; i < n; ++i)
        prod *= U(i, i);
    return prod*(luObj.perm().parity() ? -1 : 1);
}
