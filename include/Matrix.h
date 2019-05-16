//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_MATRIX_H
#define LU_MATRIX_H

#include <iterator>
#include "Vector.h"
#include "aliases.h"

typedef size_t index_t;


/// Representation of a square matrix
class Matrix {
public:
    typedef Vector Row;
    typedef std::valarray<Row> Data;
    typedef std::initializer_list<std::initializer_list<double>> InitList;

    /// Create an n*n dense matrix
    explicit Matrix(size_t n) : _n(n), _data(Row(0.0, _n), _n) {}
    Matrix(double **mat, size_t n); ///< Copy from array of pointers to rows
    Matrix(const InitList &init); ///< Construct a matrix from a 2D init list
    Matrix() = default;

    /// Assign contents from init list:
    Matrix &operator=(const InitList &init);

    // Getters:
    inline size_t size() const { return _n; }
    inline const Data &data() { return _data; }

    // Subscript operators (const and non-const). operator() has bounds checking.
    double operator()(index_t i, index_t j) const;
    double &operator()(index_t i, index_t j);
    const Row &operator[](index_t i) const { return _data[i]; }
    Row &operator[](index_t i) { return _data[i]; }

    // Operations
    Matrix &operator+=(const Matrix &other);
    Matrix &operator-=(const Matrix &other);

    // Iterator stuff:
    class iterator;
    class const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

private:
    size_t _n = 0;  ///< dimension of matrix
    Data _data;  ///< array of `std::valarray`s containing the matrix data

    // Throws exception if (i, j) is out-of-bounds (i.e. i >= _n or j >= _n)
    void checkMatrixBounds(index_t i, index_t j) const;


    class base_iterator;
};


//-------------------- Matrix operations -------------------//

Vector operator*(const Matrix &A, const Vector &v);
Matrix operator*(const Matrix &A, const Matrix &B);
Matrix operator+(const Matrix &A, const Matrix &B);
Matrix operator-(const Matrix &A, const Matrix &B);


//-------------------- Matrix::iterator -------------------//

class Matrix::base_iterator
        : public std::iterator<std::random_access_iterator_tag, double> {
public:
    base_iterator &operator++();
    base_iterator &operator--();
    base_iterator &operator+=(ptrdiff_t offset);
    base_iterator &operator-=(ptrdiff_t offset) { return *this += -offset; }
    ptrdiff_t operator-(const base_iterator &rhs) const;

    bool operator==(const base_iterator &rhs) const;
    bool operator!=(const base_iterator &rhs) const { return not(*this == rhs); }
    bool operator<(const base_iterator &rhs) const { return rhs - *this > 0; }
    bool operator>(const base_iterator &rhs) const { return *this - rhs > 0; }
    bool operator<=(const base_iterator &rhs) const { return not(*this > rhs); }
    bool operator>=(const base_iterator &rhs) const { return not(*this < rhs); }

protected:
    base_iterator(const Matrix *mat, index_t i, index_t j);

    const Matrix *mat;
    index_t i, j;
};


class Matrix::const_iterator : public Matrix::base_iterator {
public:
    explicit
    const_iterator(const Matrix *mat, index_t i = 0, index_t j = 0) : base_iterator(mat, i, j) {}

    double operator*() const { return mat->operator()(i, j); };
};


class Matrix::iterator : public Matrix::base_iterator {
public:
    explicit
    iterator(Matrix *mat, index_t i = 0, index_t j = 0) : base_iterator(mat, i, j) {}

    double &operator*() { return const_cast<Matrix *>(mat)->operator()(i, j); };
};


template<typename MatrixIter>
MatrixIter operator+(const MatrixIter &lhs, ptrdiff_t rhs) {
    MatrixIter temp = lhs;
    temp += rhs;
    return temp;
}

template<typename MatrixIter>
MatrixIter operator-(const MatrixIter &lhs, ptrdiff_t rhs) {
    MatrixIter temp = lhs;
    temp -= rhs;
    return temp;
}


#endif //LU_MATRIX_H
