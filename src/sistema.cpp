#include <cassert>

#include "resol.h"
#include "sistema.h"
#include "errors.h"
#include "norm.h"


SolveResult solve(const Matrix &A, const Vector &b, double tol) {
    try {
        LUDecomposition luObj(A, tol);
        auto &&res = SolveResult(true, std::move(luObj), solveLU(luObj, b));
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

double residue(const Matrix &A, const Vector &x, const Vector &b) {
    return normInf(A*x - b)/normInf(x);
}

ExtraSolveInfo getExtraSolveInfo(const Matrix &A, const Vector &b, const SolveResult &res) {
    const auto& luObj = res.getLU();
    return {
            .residue = residue(A, res.solution(), b),
            .cond1 = conditionNumber(A, luObj, NormType::L1),
            .condInf = conditionNumber(A, luObj, NormType::Inf)
    };
}


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
