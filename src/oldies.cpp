//
// Created by Paolo on 14/05/2019.
//

#include "oldies.h"

#include "Matrix.h"
#include "LUDecomposition.h"
#include "errors.h"

void copy(const Matrix &mat, double **a) {
    const size_t n = mat.size();
    for (index_t i = 0; i < n; ++i)
        for (index_t j = 0; j < n; ++j)
            a[i][j] = mat[i][j];
}

int lu(double **a, int n, int perm[], double tol) {
    try {
        LUDecomposition luObj(Matrix{a, size_t(n)}, tol);
        const auto &luPerm = luObj.perm();
        copy(luObj.decompMatrix(), a);
        std::copy(begin(luPerm.vector()), end(luPerm.vector()), perm);
        return (luPerm.parity() ? -1 : 1);
    } catch (SingularMatrixError &) {
        freemat(a, n);
        return 0;
    }
}

