#ifndef LU_SISTEMA_CPP
#define LU_SISTEMA_CPP

#define RESOL_DECLARATIONS_ONLY
#define LU_DECLARATIONS_ONLY

#include "resol.cpp"
#include "lu.cpp"


//---------- DECLARATIONS ----------//

struct SolveResult {
    bool success;
    Vector solution;

    inline explicit
    operator bool() const { return success; }
};

SolveResult solve(const Matrix &A, const Vector &b, double tol = numcomp::DEFAULT_TOL);


//----- C-style interface -----//

int sistema(double **a, double x[], double b[], int n, double tol);



#ifndef SISTEMA_DECLARATIONS_ONLY
//---------- IMPLEMENTATION ----------//

SolveResult solve(const Matrix &A, const Vector &b, double tol) {
    try {
        LUDecomposition luObj(A, tol);
        return {true, solve(luObj, b)};
    } catch (SingularMatrixError &) { return {false, {}}; }
}


//----- C-style interface -----//

int sistema(double **a, double *x, double *b, int n, double tol) {
    Matrix mat(a, n);
    Vector b_(b, n);
    SolveResult result = solve(mat, b_, tol);
    if (not result) return 0;
    std::copy(begin(result.solution), end(result.solution), x);
}




#endif  // #ifndef SISTEMA_DECLARATIONS_ONLY

#endif  // #ifndef LU_SISTEMA_CPP
