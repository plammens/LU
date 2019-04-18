#include <utility>

#ifndef LU_SISTEMA_CPP
#define LU_SISTEMA_CPP

#define RESOL_DECLARATIONS_ONLY
#define LU_DECLARATIONS_ONLY

#include "resol.cpp"
#include "lu.cpp"


//---------- DECLARATIONS ----------//

class SolveResult {
public:
    inline explicit
    operator bool() const { return _success; }

    const Vector &solution() const;
    inline const LUDecomposition &luDecomp() const { return _luObj; }

private:
    bool _success = false;
    LUDecomposition _luObj = {};
    Vector _solution = {};

    explicit SolveResult(bool success, Vector &&solution, LUDecomposition &&luObj);
    SolveResult() = default;

    friend SolveResult solve(const Matrix &, const Vector &, double tol);
};


SolveResult solve(const Matrix &A, const Vector &b, double tol = numcomp::DEFAULT_TOL);


//----- C-style interface -----//

int sistema(double **a, double x[], double b[], int n, double tol);




#ifndef SISTEMA_DECLARATIONS_ONLY
//---------- IMPLEMENTATION ----------//

SolveResult solve(const Matrix &A, const Vector &b, double tol) {
    try {
        LUDecomposition luObj(A, tol);
        return SolveResult(true, solve(luObj, b), std::move(luObj));
    } catch (SingularMatrixError &) {
        return SolveResult();
    }
}


SolveResult::SolveResult(bool success, Vector &&solution, LUDecomposition &&luObj)
        : _success(success), _luObj(std::move(luObj)), _solution(std::move(solution)) {}

const Vector &SolveResult::solution() const {
    if (not _success)
        throw std::invalid_argument("invalid request for unsuccessful solve's solution");
    return _solution;
}


//----- C-style interface -----//

int sistema(double **a, double *x, double *b, int n, double tol) {
    Matrix mat(a, n);
    Vector b_(b, n);
    SolveResult result = solve(mat, b_, tol);
    if (not result) return 0;
    const Vector &solution = result.solution();
    std::copy(begin(solution), end(solution), x);
    return result.luDecomp().perm().parity()? -1 : 1;
}




#endif  // #ifndef SISTEMA_DECLARATIONS_ONLY

#endif  // #ifndef LU_SISTEMA_CPP
