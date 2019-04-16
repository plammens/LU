#ifndef LU_LU_CPP
#define LU_LU_CPP

#include <valarray>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <cmath>



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

template<typename T>
using nested_init_list = std::initializer_list<std::initializer_list<T>>;


/// Representation of a square matrix
class Matrix {
public:
    /// Create an n*n dense matrix with entries set to 0
    explicit Matrix(size_t n) : Matrix(n, n*n) {}
    /// Construct a matrix from a 2-dimensional init list
    Matrix(const nested_init_list<double> &init);

    // Getters:
    inline size_t size() const { return _n; }
    inline const std::valarray<double> &data() { return _data; }

    // Subscript operators (const and non-const):
    virtual double operator()(index_t i, index_t j) const;
    virtual double &operator()(index_t i, index_t j);

    // Iterator stuff:
    class iterator;

    class const_iterator;

    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();

protected:
    size_t _n;  /// dimension of matrix
    std::valarray<double> _data;  /// flat array containing the data

    // Constructs an n*n Matrix with custom data size (intended for subclasses)
    Matrix(size_t dimension, size_t dataSize) : _n(dimension), _data(0.0, dataSize) {}

    // Throws exception if (i, j) is out-of-bounds (i.e. i >= _n or j >= _n)
    void checkMatrixBounds(index_t i, index_t j) const;

    // Get the appropriate index for the internal flat array
    virtual inline
    index_t flattenIndex(index_t i, index_t j) const { return _n*i + j; };

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
    base_iterator(const Matrix *mat, index_t i, index_t j) : mat(mat), i(i), j(j) {}

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



//----- Triangular -----//

/// Abstract Matrix subclass for triangular matrices
class Triangular : public Matrix {
public:
    double operator()(index_t i, index_t j) const override;
    double &operator()(index_t i, index_t j) override;

protected:
    /// Construct the internals of an n by n triangular matrix
    explicit Triangular(size_t n) : Matrix(n, n*(n + 1)/2) {}

    /// Initialize the contents from an 2D init list
    void initFromList(const nested_init_list<double> &init);

    /// Whether the element at (i, j) is inside the "triangle"
    /// (defined by the particular subclass)
    virtual bool isInTriangle(index_t i, index_t j) const = 0;

    /// Size of the ith row of the triangle
    virtual size_t rowSize(index_t i) const = 0;
};


class LowerTriangular : public Triangular {
public:
    explicit LowerTriangular(size_t n) : Triangular(n) {}
    LowerTriangular(const nested_init_list<double> &init) : Triangular(init.size()) {
        initFromList(init);
    }

private:
    inline index_t flattenIndex(index_t i, index_t j) const override { return i*(i + 1)/2 + j; }
    inline bool isInTriangle(index_t i, index_t j) const override { return j <= i; }
    inline size_t rowSize(index_t i) const override { return i; }
};


class UpperTriangular : public Triangular {
public:
    explicit UpperTriangular(size_t n) : Triangular(n) {}
    UpperTriangular(const nested_init_list<double> &init) : Triangular(init.size()) {
        initFromList(init);
    }

private:
    inline index_t flattenIndex(index_t i, index_t j) const override { return i*_n - i*(i + 1)/2 + j; }
    inline bool isInTriangle(index_t i, index_t j) const override { return j >= i; }
    inline size_t rowSize(index_t i) const override { return _n - i; }
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
    LUDecomposition(Matrix &mat, double tol = numcomp::DEFAULT_TOL);

    /// Get the permutation vector
    const std::valarray<index_t> &perm() const { return _perm; }

private:
    Matrix _decomp_mat;  ///< decomposition matrix (internal data storage)
    std::valarray<index_t> _perm;  ///< permutation vector
    double _tol;  ///< numerical tolerance

    /// Performs the actual LU decomposition. Called upon construction.
    void decompose();
    /// Swaps the row at pivot_index with the largest (relatively scaled) pivot
    void scaledPartialPivoting(index_t pivot_index);
};


#ifndef LU_DECLARATIONS_ONLY

//--------------- IMPLEMENTATION ---------------//

//----- Matrix -----//

Matrix::Matrix(const nested_init_list<double> &init) : _n(init.size()), _data(_n*_n) {
    index_t i = 0;
    for (const auto &row : init) {
        if (row.size() != _n)
            throw std::invalid_argument("unevenly sized init list for square matrix");
        for (double num : row) _data[i++] = num;
    }
}

double Matrix::operator()(index_t i, index_t j) const {
    checkMatrixBounds(i, j);
    return _data[flattenIndex(i, j)];
}

double &Matrix::operator()(index_t i, index_t j) {
    checkMatrixBounds(i, j);
    return _data[flattenIndex(i, j)];
}

inline
void Matrix::checkMatrixBounds(index_t i, index_t j) const {
    if (i >= _n or j >= _n) throw std::out_of_range("matrix subscript out of range");
}

Matrix::const_iterator Matrix::begin() const {
    return const_iterator(this, 0, 0);
}

Matrix::const_iterator Matrix::end() const {
    return const_iterator(this, _n, 0);
}

Matrix::iterator Matrix::begin() {
    return iterator(this, 0, 0);
}

Matrix::iterator Matrix::end() {
    return iterator(this, _n, 0);
}

//----- Matrix::iterator -----//

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
    return mat->flattenIndex(i, j) - mat->flattenIndex(rhs.i, rhs.j);
}

bool Matrix::base_iterator::operator==(const Matrix::base_iterator &rhs) const {
    return mat == rhs.mat and i == rhs.i and j == rhs.j;
}


//----- Triangular -----//

double Triangular::operator()(index_t i, index_t j) const {
    checkMatrixBounds(i, j);
    if (not isInTriangle(i, j)) return 0;
    return _data[flattenIndex(i, j)];
}

double &Triangular::operator()(index_t i, index_t j) {
    checkMatrixBounds(i, j);
    if (not isInTriangle(i, j))
        throw std::out_of_range("cannot write to null side of triangular matrix");
    return _data[flattenIndex(i, j)];
}

void Triangular::initFromList(const nested_init_list<double> &init) {
    index_t flattened = 0;

    index_t i = 0;
    for (const auto &initRow : init) {
        size_t initRowSize = initRow.size();
        if (initRowSize < rowSize(i) or initRowSize > _n)
            throw std::invalid_argument("unevenly sized init list");

        index_t j = 0;
        for (double num : initRow) {
            if (not isInTriangle(i, j)) {
                if (num != 0.0)
                    throw std::invalid_argument("null side of triangular matrix should be zero");
            } else _data[flattened++] = num;
            ++j;
        }

        ++i;
    }
}


//---------- functions ----------//

/**
 * Returns the maximal absolute value of the numbers in a range.
 * @param first,last pair of forward input iterators describing the range
 */
template<class FIter>
double max_abs(FIter first, FIter last) {
    double max_value = 0;
    for (; first != last; ++first)
        max_value = std::max(std::abs(*first), max_value);
    return max_value;
}



//----- LUDecomposition -----//

LUDecomposition::LUDecomposition(const Matrix &mat, double tol)
        : _decomp_mat(mat), _tol(tol), _perm(mat.size()) {
    decompose();
}

void LUDecomposition::decompose() {
    const size_t n = mat.size();
    for (index_t pivot_index = 0; pivot_index < n - 1; ++pivot_index) {
        // Swap rows if necessary, and get pivot:
        scaledPartialPivoting(pivot_index);
        const double pivot = _decomp_mat(pivot_index, pivot_index);

        // Update rows below:
        for (index_t i = pivot_index + 1; i < n; ++i) {
            const double multiplier = -mat(i, pivot_index)/pivot;
            _decomp_mat(i, pivot_index) = 0.0;  // shortcut
            for (index_t j = pivot_index + 1; j < n; ++j) // add pivot row to ith row
                _decomp_mat(i, j) += multiplier*mat(pivot_index, j);
            _decomp_mat(i, pivot_index) = multiplier;  // write multiplier to lower triangle
        }
    }
}

void LUDecomposition::scaledPartialPivoting(index_t pivot_index) {
    size_t n = _decomp_mat.size();
    struct { double value; Index index; } max_pivot = {0, pivot_index};

    for (Index i = pivot_index; i < n; ++i) {
        const Matrix::Row &row = mat[i];
        double max_abs_value = max_abs(row.begin() + pivot_index, row.end());  // scaling factor
        if (numeric::isnull(max_abs_value, tol)) return {false};  // null row --> singular matrix
        // update max pivot:
        double scaled_pivot = row[pivot_index]/max_abs_value;
        if (scaled_pivot > max_pivot.value) max_pivot = {scaled_pivot, i};
    }

    std::swap(mat[pivot_index], mat[max_pivot.index]);  // constant complexity; just swaps pointers
    return {true, perm};
}

//int lu(double **mat, int n, int perm[], double tol) {
//}


#endif  // LU_DECLARATIONS_ONLY

#endif  // LU_LU_CPP
