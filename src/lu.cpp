#ifndef LU_LU_CPP
#define LU_LU_CPP

#include <valarray>
#include <initializer_list>
#include <iterator>
#include <exception>
#include <cmath>
#include <cassert>



//--------------- DECLARATIONS ---------------//

/// Namespace for numerical comparisons
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


// Aliases and using declarations:
using std::size_t;
using std::ptrdiff_t;
typedef size_t index_t;
using std::begin;
using std::end;

/// Representation of a square matrix
class Matrix {
public:
    typedef std::valarray<double> Row;
    typedef std::valarray<Row> Data;
    typedef std::initializer_list<std::initializer_list<double>> InitList;

    /// Create an n*n dense matrix
    explicit Matrix(size_t n) : _n(n), _data(Row(_n), _n) {}

    /// Construct a matrix from a 2-dimensional init list
    Matrix(const InitList &init);
    Matrix &operator=(const InitList &init);

    // Getters:
    inline size_t size() const { return _n; }
    inline const Data &data() { return _data; }

    // Subscript operators (const and non-const):
    virtual double operator()(index_t i, index_t j) const;
    virtual double &operator()(index_t i, index_t j);
    const Row &operator[](index_t i) const { return _data[i]; }
    Row &operator[](index_t i) { return _data[i]; }

    // Iterator stuff:
    class iterator;
    class const_iterator;
    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

protected:
    size_t _n;  /// dimension of matrix
    Data _data;  /// flat array containing the data

    // Throws exception if (i, j) is out-of-bounds (i.e. i >= _n or j >= _n)
    void checkMatrixBounds(index_t i, index_t j) const;

private:
    class base_iterator;
};

//----- Matrix::iterator -----//

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




//----- LU -----//

class LUDecomposition {
public:
    /**
     * Compute the LU decomposition of a matrix.
     * @param mat  matrix to decompose
     * @param tol  numerical tolerance
     * @throws SingularMatrixError if mat is singular
     */
    explicit
    LUDecomposition(const Matrix &mat, double tol = numcomp::DEFAULT_TOL);

    // getters:
    const std::valarray<index_t> &perm() const { return _perm; }
    const Matrix &getDecompMatrix() { return _mat; }

private:
    Matrix _mat;  ///< decomposition matrix (internal data storage)
    std::valarray<index_t> _perm;  ///< permutation vector
    double _tol;  ///< numerical tolerance

    /// Performs the actual LU decomposition. Called upon construction.
    void decompose();

    /// Swaps the row at pivot_index with the row with the (relatively) largest pivot.
    void scaledPartialPivoting(index_t pivot_index);
};


/// Indicates that an algorithm encountered a singular matrix
class SingularMatrixError : public std::exception {
    static constexpr auto message = "singular matrix";
public:
    const char *what() const noexcept override { return message; }
};





#ifndef LU_DECLARATIONS_ONLY

//--------------- IMPLEMENTATION ---------------//

//----- Matrix -----//

Matrix::Matrix(const Matrix::InitList &init) : Matrix(init.size()) {
    this->operator=(init);
}

Matrix& Matrix::operator=(const InitList &init) {
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
    if (++j == mat->size()) { ++i; j = 0; }
    return *this;
}

Matrix::base_iterator &Matrix::base_iterator::operator+=(ptrdiff_t offset) {
    size_t n = mat->size();
    i += (j += offset)/n;
    j %= n;
    return *this;
}

Matrix::base_iterator &Matrix::base_iterator::operator--() {
    if (j == 0) { --i; j = mat->size() - 1; }
    else --j;
    return *this;
}

ptrdiff_t Matrix::base_iterator::operator-(const Matrix::base_iterator &rhs) const {
    return mat->_n*(ptrdiff_t(i) - ptrdiff_t(rhs.i)) + ptrdiff_t(j) - ptrdiff_t(rhs.j);
}

bool Matrix::base_iterator::operator==(const Matrix::base_iterator &rhs) const {
    return mat == rhs.mat and i == rhs.i and j == rhs.j;
}




//----- LUDecomposition -----//

LUDecomposition::LUDecomposition(const Matrix &mat, double tol)
        : _mat(mat), _perm(mat.size()), _tol(tol) {
    // initialize permutation vector:
    for (index_t i = 0; i < mat.size(); ++i) _perm[i] = i;
    decompose();
}

void LUDecomposition::decompose() {
    const size_t n = _mat.size();

    for (index_t pivot_index = 0; pivot_index < n - 1; ++pivot_index) {
        // Swap rows if necessary, and get pivot:
        scaledPartialPivoting(pivot_index);
        const Matrix::Row &pivot_row = _mat[pivot_index];
        const double pivot = pivot_row[pivot_index];

        // Update rows below the pivot row:
        for (index_t i = pivot_index + 1; i < n;++i) {
            Matrix::Row &row = _mat[i];
            const double multiplier = row[pivot_index]/pivot;
            // Subtract multiple of pivot row from current row:
            auto slice = std::slice(pivot_index + 1, n - pivot_index - 1, 1);
            row[slice] -= multiplier*pivot_row[slice];
            row[pivot_index] = multiplier;  // write multiplier to lower triangle
        }
    }
}

/* Helper function.
 * Returns the maximal absolute value of the numbers in a range.
 * @param first,last pair of forward input iterators describing the range
 * @return the maximal absolute value
 */
template<class FIter>
double max_abs(FIter first, FIter last) {
    double max_value = 0;
    for (; first != last; ++first)
        max_value = std::max(std::abs(*first), max_value);
    return max_value;
}

void LUDecomposition::scaledPartialPivoting(index_t pivot_index) {
    const size_t n = _mat.size();
    struct { double value; index_t index; } max_pivot = {0, pivot_index};

    for (index_t i = pivot_index; i < n; ++i) {
        const auto &row = _mat[i];
        double max_abs_value = max_abs(begin(row) + pivot_index, end(row));  // scaling factor
        if (numcomp::isnull(max_abs_value, _tol)) throw SingularMatrixError();  // null row --> singular matrix
        // update max pivot:
        double scaled_pivot = row[pivot_index]/max_abs_value;
        if (scaled_pivot > max_pivot.value) max_pivot = {scaled_pivot, i};
    }

    // Swap the indices with the row with maximal pivot:
    std::swap(_perm[pivot_index], _perm[max_pivot.index]);
    std::swap(_mat[pivot_index], _mat[max_pivot.index]);  // constant complexity; just swaps pointers
}

//int lu(double **mat, int n, int perm[], double tol) {
//}


#endif  // #ifndef LU_DECLARATIONS_ONLY

#endif  // #ifndef LU_LU_CPP
