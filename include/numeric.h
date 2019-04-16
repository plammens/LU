/// @file
/// Inline functions for numerical comparisons (header-only)

#ifndef LU_NUMERIC_H
#define LU_NUMERIC_H

#include <cmath>

namespace numcomp {
    
    constexpr const double DEFAULT_TOL = 1e-12;
    
    inline
    bool equal(double a, double b, double tol = DEFAULT_TOL) {
        return std::abs(a - b) < tol;
    }

    inline
    bool isnull(double a, double tol = DEFAULT_TOL) {
        return equal(a, 0, tol);
    }

    inline
    bool leq(double a, double b, double tol = DEFAULT_TOL) {
        return a < b + tol;
    }

    inline
    bool geq(double a, double b, double tol = DEFAULT_TOL) {
        return leq(b, a, tol);
    }

    inline
    bool less(double a, double b, double tol = DEFAULT_TOL) {
        return not leq(b, a, tol);
    }

    inline
    bool greater(double a, double b, double tol = DEFAULT_TOL) {
        return less(b, a, tol);
    }

}

#endif //LU_NUMERIC_H
