//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_RESOL_H
#define LU_RESOL_H

#include "Vector.h"
#include "LUDecomposition.h"

/// Compute the solution of a linear system given the LU decomposition of the matrix
Vector solveLU(const LUDecomposition &luObj, const Vector &b);

/// C-style interface to solve(const LUDecomposition &, const Vector &)
void resol(double **a, double x[], double b[], int n, int perm[]);

#endif //LU_RESOL_H
