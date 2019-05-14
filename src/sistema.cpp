#include <cassert>

#include "resol.h"
#include "sistema.h"
#include "errors.h"


SolveResult solve(const Matrix &A, const Vector &b, double tol) {
    try {
        LUDecomposition luObj(A, tol);
        auto &&res = SolveResult(true, std::move(luObj), solveLU(luObj, b));
        res._residue = residue(A, res._solution, b);
        return std::move(res);
    } catch (SingularMatrixError &) {
        return SolveResult(tol);
    }
}


SolveResult::SolveResult(bool success, LUDecomposition &&luObj, Vector &&solution)
        : _success(success), _luObj(std::move(luObj)), _solution(std::move(solution)) {}

const Vector &SolveResult::solution() const {
    if (not _success)
        throw std::invalid_argument("invalid request for unsuccessful solve's solution");
    return _solution;
}

Vector operator*(const Matrix &A, const Vector &v) {
    const size_t n = A.size();
    assert(n == v.size());

    Vector result(0.0, n);
    for (index_t i = 0; i < n; ++i) {
        double &elem = result[i];
        for (index_t j = 0; j < n; ++j)
            elem += A[i][j]*v[j];
    }

    return result;
}

inline
double norm(const Vector &v) {
    return std::abs(v).max();
}

double residue(const Matrix &A, const Vector &x, const Vector &b) {
    return norm(A*x - b)/norm(x);
}


//double conditionNumber(const Matrix &A, NormType nt) {
//    // TODO: implement inverse
//    switch (nt) {
//        case NormType::L1: return norm1(A);
//        case NormType::Inf: return normInf(A);
//    }
//}


//double norm1(const Matrix& A) {
//    size_t n = A.size();
//    Vector colSums(0.0, n);
//
//    for (const Matrix::Row &row : A)
//        colSums += std::abs(row);
//
//    return colSums.max();
//}
//
//
//double normInf(const Matrix&A) {
//    size_t n = A.size();
//    double maxSum = 0.0;
//
//    for (const Matrix::Row &row : A)
//        maxSum = std::max(maxSum, std::abs(row).sum());
//
//    return maxSum;
//}



//----- C-style interface -----//

int sistema(double **a, double *x, double *b, int n, double tol) {
    Matrix mat(a, n);
    Vector b_(b, n);
    SolveResult result = solve(mat, b_, tol);
    if (not result) return 0;
    const Vector &solution = result.solution();
    std::copy(begin(solution), end(solution), x);
    return result.getLU().perm().parity() ? -1 : 1;
}
