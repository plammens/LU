//
// Created by Paolo on 14/05/2019.
//

#include "Matrix.h"


Matrix::Matrix(const Matrix::InitList &init) : Matrix(init.size()) {
    this->operator=(init);
}

Matrix::Matrix(double **mat, size_t n) : _n(n), _data(_n) {
    for (index_t i = 0; i < n; ++i)
        _data[i] = std::valarray<double>(mat[i], n);
}

Matrix &Matrix::operator=(const InitList &init) {
    index_t i = 0;
    for (const auto &initRow : init) {
        if (initRow.size() != _n)
            throw std::invalid_argument("badly shaped init list for square matrix");
        auto &row = _data[i];
        index_t j = 0;
        for (double num : initRow) row[j++] = num;
        ++i;
    }
    return *this;
}

double Matrix::operator()(index_t i, index_t j) const {
    checkMatrixBounds(i, j);
    return _data[i][j];
}

double &Matrix::operator()(index_t i, index_t j) {
    checkMatrixBounds(i, j);
    return _data[i][j];
}

inline
void Matrix::checkMatrixBounds(index_t i, index_t j) const {
    if (i >= _n or j >= _n) throw std::out_of_range("matrix subscript out of range");
}

Matrix::const_iterator Matrix::begin() const { return const_iterator(this, 0, 0); }
Matrix::const_iterator Matrix::end() const { return const_iterator(this, _n, 0); }
Matrix::iterator Matrix::begin() { return iterator(this, 0, 0); }
Matrix::iterator Matrix::end() { return iterator(this, _n, 0); }



//----- Matrix::iterator -----//

Matrix::base_iterator::base_iterator(const Matrix *mat, index_t i, index_t j)
        : mat(mat), i(i), j(j) {
    const size_t n = mat->size();
    this->i += j/n;
    this->j %= n;
}

Matrix::base_iterator &Matrix::base_iterator::operator++() {
    if (++j == mat->size()) {
        ++i;
        j = 0;
    }
    return *this;
}

Matrix::base_iterator &Matrix::base_iterator::operator+=(ptrdiff_t offset) {
    size_t n = mat->size();
    i += (j += offset)/n;
    j %= n;
    return *this;
}

Matrix::base_iterator &Matrix::base_iterator::operator--() {
    if (j == 0) {
        --i;
        j = mat->size() - 1;
    }
    else --j;
    return *this;
}

ptrdiff_t Matrix::base_iterator::operator-(const Matrix::base_iterator &rhs) const {
    return mat->_n*(ptrdiff_t(i) - ptrdiff_t(rhs.i)) + ptrdiff_t(j) - ptrdiff_t(rhs.j);
}

bool Matrix::base_iterator::operator==(const Matrix::base_iterator &rhs) const {
    return mat == rhs.mat and i == rhs.i and j == rhs.j;
}

