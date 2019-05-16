//
// Created by Paolo on 16/05/2019.
//

#include "norm.h"

#include <cmath>
#include <functional>
#include "errors.h"
#include "inverse.h"


VectorNorm getVectorNorm(NormType nt) {
    switch (nt) {
        case NormType::L1: return norm1;
        case NormType::L2: return norm2;
        case NormType::Inf: return normInf;
    }
}

MatrixNorm getMatrixNorm(NormType nt) {
    switch (nt) {
        case NormType::L1: return norm1;
        case NormType::Inf: return normInf;
        default: throw ValueError("unknown matrix norm");
    }
}


double norm(const Vector &vec, NormType nt) {
    return getVectorNorm(nt)(vec);
}

double norm1(const Vector &vec) {
    return std::abs(vec).sum();
}

double norm2(const Vector &vec) {
    return std::sqrt((vec*vec).sum());
}

double normInf(const Vector &vec) {
    return std::abs(vec).max();
}

double norm(const Matrix &A, NormType nt) {
    return getMatrixNorm(nt)(A);
}

double norm1(const Matrix &A) {
    // Max column absolute sum
    const size_t n = A.size();
    Vector colSums(0.0, n);
    for (index_t i = 0; i < n; ++i) colSums += std::abs(A[i]);
    return colSums.max();
}

double normInf(const Matrix &A) {
    // Max row absolute sum
    const size_t n = A.size();
    double maxSum = 0.0;
    for (index_t i = 0; i < n; ++i)
        maxSum = std::max(maxSum, std::abs(A[i]).sum());
    return maxSum;
}

double conditionNumber(const Matrix &A, const Matrix &AInv, NormType nt) {
    auto matNorm = getMatrixNorm(nt);
    return matNorm(A)*matNorm(AInv);
}

double conditionNumber(const Matrix &A, const LUDecomposition &luObj, NormType nt) {
    // TODO: inverse
    return conditionNumber(A, inverse(luObj), nt);
}

double conditionNumber(const Matrix &A, NormType nt) {
    return conditionNumber(A, LUDecomposition(A), nt);
}
