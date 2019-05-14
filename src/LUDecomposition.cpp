//
// Created by Paolo on 14/05/2019.
//

#include "LUDecomposition.h"
#include "errors.h"

//----- LUDecomposition -----//

LUDecomposition::LUDecomposition(const Matrix &mat, double tol)
        : _mat(mat), _perm(mat.size()), _tol(tol) {
    decompose();
}

void LUDecomposition::decompose() {
    const size_t n = _mat.size();

    for (index_t pivot_index = 0; pivot_index < n; ++pivot_index) {
        // Swap rows if necessary, and get pivot:
        scaledPartialPivoting(pivot_index);
        const Matrix::Row &pivot_row = _mat[pivot_index];
        const double pivot = pivot_row[pivot_index];

        // Update rows below the pivot row:
        for (index_t i = pivot_index + 1; i < n; ++i) {
            Matrix::Row &row = _mat[i];
            double &multiplier = row[pivot_index] /= pivot;
            // Subtract multiple of pivot row from current row:
            auto slice = std::slice(pivot_index + 1, n - pivot_index - 1, 1);
            row[slice] -= multiplier*pivot_row[slice];
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
        // if max_abs_value is zero, the whole row is, so the matrix is singular:
        if (max_abs_value < _tol) throw SingularMatrixError(_tol);
        // update max pivot:
        double scaled_pivot = std::abs(row[pivot_index])/max_abs_value;
        if (scaled_pivot > max_pivot.value) max_pivot = {scaled_pivot, i};
    }

    // Swap the indices with the row with maximal pivot:
    _perm.permute(pivot_index, max_pivot.index);
    std::swap(_mat[pivot_index], _mat[max_pivot.index]);  // constant complexity; just swaps pointers
}

