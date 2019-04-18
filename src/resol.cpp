#ifndef LU_RESOL_CPP
#define LU_RESOL_CPP

#define LU_DECLARATIONS_ONLY
#include "lu.cpp"

//---------- DECLARATIONS ----------//

typedef std::valarray<double> Vector;

Vector solve(const LUDecomposition &luObj, const Vector &b);


//----- C-style interface -----//

void resol(double ** a , double x[] , double b[] , int n , int perm[]);




#ifndef RESOL_DECLARATIONS_ONLY
//---------- IMPLEMENTATION ----------//

Vector solveLower(const Matrix &lower, const Vector &b, const Permutation::Vector &perm) {
    const size_t n = lower.size();

    Vector y(n);
    for (index_t i = 0; i < n; ++i) {
        double &elem = y[i] = b[perm[i]];
        for (index_t j = 0; j < i; ++j) elem -= lower[i][j]*y[j];
    }

    return y;
}

Vector solveUpper(const Matrix &upper, const Vector &y) {
    const size_t n = upper.size();

    Vector x(n);
    for (long i = n - 1; i >= 0; --i) {
        double &elem = x[i] = y[i];
        for (index_t j = i + 1; j < n; ++j) elem -= upper[i][j]*x[j];
        elem /= upper[i][i];
    }

    return x;
}


Vector solve(const LUDecomposition &luObj, const Vector &b) {
    const Matrix &decompMat = luObj.decompMatrix();
    const auto &perm = luObj.perm().vector();
    if (decompMat.size() != b.size())
        throw std::invalid_argument("vector and matrix size are different");

    const Vector &&y = solveLower(decompMat, b, perm);
    return solveUpper(decompMat, y);
}



//----- C-style interface -----//

void resol(double **a, double *x, double *b, int n, int *perm) {
    Matrix mat(a, n);
    Permutation::Vector perm_(n);
    std::copy(perm, perm + n, begin(perm_));
    Vector b_(b, n);

    const Vector &&y_ = solveLower(mat, b_, perm_);
    const Vector &&x_ = solveUpper(mat, y_);
    std::copy(begin(x_), end(x_), x);
}


#endif  // #ifndef RESOL_DECLARATIONS_ONLY

#endif // #ifndef LU_RESOL_CPP

