#include <utility>

#ifndef LU_SISTEMA_CPP
#define LU_SISTEMA_CPP

#define RESOL_DECLARATIONS_ONLY
#define LU_DECLARATIONS_ONLY

#include "resol.cpp"
#include "lu.cpp"


//---------- DECLARATIONS ----------//

/// Utility class to store the result of solving a linear system
class SolveResult {
public:
    /// Returns whether solving was successful
    inline explicit
    operator bool() const { return _success; }

    /**
     * Numerical solution vector.
     * @return the computed solution to the linear system
     * @pre `bool(this)` is `true`.
    */
    const Vector &solution() const;

    inline const LUDecomposition &getLU() const { return _luObj; }
    inline double tol() const { return _tol; }
    inline double residue() const { return _residue; }

private:
    bool _success = false;
    LUDecomposition _luObj = {};
    Vector _solution = {};
    double _tol = 0.0;
    double _residue = 0.0;

    explicit SolveResult(bool success, Vector &&solution, LUDecomposition &&luObj);
    SolveResult() = default;

    friend SolveResult solve(const Matrix &, const Vector &, double tol);
};


/**
 * Solve the linear system Ax = b
 * @param A  matrix
 * @param b  vector of independent terms
 * @return  a SolveResult object which evaluates to `true` if
 * the procedure was successful, and `false` otherwise (A is singular).
 * In the former case, result.solution() will be the numerical solution to
 * the system Ax = b.
 */
SolveResult solve(const Matrix &A, const Vector &b, double tol = numcomp::DEFAULT_TOL);

/**
 * Calculate the residue for the approximate solution x to the system Ax = b
 * @param A  matrix
 * @param b  vector of independent terms
 * @param x  approximate solution
 * @return  ||Ax - b||/||x||
 */
double residue(const Matrix &A, const Vector &b, const Vector &x);



//----- C-style interface -----//

int sistema(double **a, double x[], double b[], int n, double tol);




#ifndef SISTEMA_DECLARATIONS_ONLY
//---------- IMPLEMENTATION ----------//

SolveResult solve(const Matrix &A, const Vector &b, double tol) {
    try {
        LUDecomposition luObj(A, tol);
        auto &&res = SolveResult(true, solve(luObj, b), std::move(luObj));
        res._tol = tol;
        res._residue = residue(A, res._solution, b);
        return std::move(res);
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

double residue(const Matrix &A, const Vector &b, const Vector &x) {
    return norm(A*x - b)/norm(x);
}


//----- C-style interface -----//

int sistema(double **a, double *x, double *b, int n, double tol) {
    Matrix mat(a, n);
    Vector b_(b, n);
    SolveResult result = solve(mat, b_, tol);
    if (not result) return 0;
    const Vector &solution = result.solution();
    std::copy(begin(solution), end(solution), x);
    return result.getLU().perm().parity()? -1 : 1;
}



#endif  // #ifndef SISTEMA_DECLARATIONS_ONLY

#endif  // #ifndef LU_SISTEMA_CPP
