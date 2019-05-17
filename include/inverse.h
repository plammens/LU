//
// Created by Paolo on 16/05/2019.
//

#ifndef LU_INVERSE_H
#define LU_INVERSE_H

#include "LUDecomposition.h"

/**
 * Matrix inverse through LU decomposition.
 * @param A  matrix to invert
 * @return  the inverse of `A` if `A` is non-singular
 * @throws  SingularMatrixError if `A` is singular
 */
Matrix inverse(const Matrix &A);

/**
 * Matrix inverse given its LU decomposition.
 * @param luObj  lu decomposition of the matrix to invert
 * @return  the inverse of the matrix
 */
Matrix inverse(const LUDecomposition &luObj);

/// n by n identity matrix
Matrix identity(size_t n);

#endif //LU_INVERSE_H
