/// @file C-style array interface
// Created by Paolo on 14/05/2019.
//

#ifndef LU_OLDIES_H
#define LU_OLDIES_H

#include "Matrix.h"
#include "aliases.h"

/// Allocate space for an n by n matrix
inline
double **newmat(size_t n) {
    auto a = new double *[n];
    for (index_t i = 0; i < n; ++i) a[i] = new double[n]();
    return a;
}

/// Free memory occupied by an n by n matrix
inline
void freemat(double **a, size_t n) {
    for (index_t i = 0; i < n; ++i) delete[] a[i];
    delete[] a;
}


/// Copy the contents of a matrix to a c-style array
void copy(const Matrix &mat, double **a);

#endif //LU_OLDIES_H
