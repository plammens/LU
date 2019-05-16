//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_SISTEMA_H
#define LU_SISTEMA_H

#include "LUDecomposition.h"
#include "Vector.h"


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
    inline double tol() const { return _success ? _luObj.tol() : _tol; }
    inline double residue() const { return _residue; }

private:
    bool _success = false;
    LUDecomposition _luObj = {};
    Vector _solution = {};
    double _residue = 0.0;
    double _tol = 0.0;

    explicit SolveResult(bool success, LUDecomposition &&luObj, Vector &&solution);
    explicit SolveResult(double tol) : _tol(tol) {}

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
 * Calculate the relative residue for the approximate solution x
 * to the linear system Ax = b
 * @param A  matrix
 * @param x  approximate solution
 * @param b  vector of independent terms
 * @return  ||Ax - b||_∞/||x||_∞
 */
double residue(const Matrix &A, const Vector &x, const Vector &b);


//----- C-style interface -----//

int sistema(double **a, double x[], double b[], int n, double tol);


#endif //LU_SISTEMA_H
