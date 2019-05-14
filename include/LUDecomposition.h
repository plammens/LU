//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_LUDECOMPOSITION_H
#define LU_LUDECOMPOSITION_H


#include "Matrix.h"
#include "Permutation.h"
#include "numcomp.h"


/// This class handles the LU decomposition of a matrix
class LUDecomposition {
public:
    /**
     * Compute the LU decomposition of a matrix.
     * @param mat  matrix to decompose
     * @param tol  numerical tolerance
     * @throws SingularMatrixError if mat is singular
     */
    explicit LUDecomposition(const Matrix &mat, double tol = numcomp::DEFAULT_TOL);
    LUDecomposition() = default;

    // getters:
    inline const Permutation &perm() const { return _perm; }
    inline const Matrix &decompMatrix() const { return _mat; }
    inline double tol() const { return _tol; }

private:
    Matrix _mat;  ///< decomposition matrix (internal data storage)
    Permutation _perm;
    double _tol = numcomp::DEFAULT_TOL;  ///< numerical tolerance

    /// Performs the actual LU decomposition. Called upon construction.
    void decompose();

    /// Swaps the row at pivot_index with the row with the (relatively) largest pivot.
    void scaledPartialPivoting(index_t pivot_index);
};


/// C-style interface
int lu(double **a, int n, int perm[], double tol);

#endif //LU_LUDECOMPOSITION_H
