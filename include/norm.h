//
// Created by Paolo on 16/05/2019.
//

#ifndef LU_NORM_H
#define LU_NORM_H

#include "LUDecomposition.h"
#include "Vector.h"

enum class NormType {
    L1 = 1,     // L1 norm
    L2 = 2,     // L2 norm
    Inf = -1,   // L_inf norm
};

typedef double (*VectorNorm)(const Vector &);
typedef double (*MatrixNorm)(const Matrix &);

// Vector norms
double norm(const Vector &vec, NormType nt);
double norm1(const Vector &vec);
double norm2(const Vector &vec);
double normInf(const Vector &vec);

// Matrix norms
double norm(const Matrix &A, NormType nt);
double norm1(const Matrix &A);
double normInf(const Matrix &A);

/**
 * Condition number of a matrix, given its inverse
 * @param A  matrix
 * @param LU  LU decomposition of `A`
 * @param nt  norm type
 * @return  the condition number of `A` under the given norm type
 *
 * @pre  `LU` is the LU decomposition of `A`
 */
double conditionNumber(const Matrix &A, NormType nt);

/**
 * Condition number of a matrix, given its inverse
 * @param A  matrix
 * @param LU  LU decomposition of `A`
 * @param nt  norm type
 * @return  the condition number of `A` under the given norm type
 *
 * @pre  `LU` is the LU decomposition of `A`
 */
double conditionNumber(const Matrix &A, const LUDecomposition &LU, NormType nt);

/**
 * Condition number of a matrix, given its inverse
 * @param A  matrix
 * @param AInv  inverse of `A`
 * @param nt  norm type
 * @return  the condition number of `A` under the given norm type
 *
 * @pre  `AInv` is the inverse of `A`
 */
double conditionNumber(const Matrix &A, const Matrix &AInv, NormType nt);


#endif //LU_NORM_H
