//
// Created by Paolo on 17/05/2019.
//

#ifndef LU_EXTRA_H
#define LU_EXTRA_H

#include "LUDecomposition.h"

/**
 * Matrix determinant through LU decomposition.
 * @param A  matrix
 * @return  the determinant of `A` up to machine precision if
 *          `A` is non-singular, otherwise 0.0
 */
double determinant(const Matrix &A);

/**
 * Matrix determinant given its LU decomposition.
 * @param luObj  lu decomposition of the matrix to invert
 * @return  the determinant of the matrix whose LU decomposition is luObj
 */
double determinant(const LUDecomposition & luObj);


#endif //LU_EXTRA_H
