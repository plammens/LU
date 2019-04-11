#include <valarray>
#include <algorithm>
#include <numeric>
#include "numeric.h"

using std::size_t;
typedef std::valarray<double> Row;
typedef std::valarray<Row> Matrix;
typedef size_t Index;
using std::begin;
using std::end;


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


/**
 * Performs scaled partial pivoting on a matrix at the row with index `pivot_index`.
 * Swaps the row at the index with the row below it with the largest (in relative absolute
 * value) pivot element.
 * Subroutine of lu(double**, int perm[], double tol).
 *
 * @param mat  matrix to perform SPP on
 * @param pivot_index  index of pivot row/column
 * @param tol  numerical tolerance
 * @return  1 on success, 0 on failure (matrix is singular)
 */
int scaled_partial_pivoting(Matrix &mat, Index pivot_index, double tol) {
    size_t n = mat.size();
    struct { double value; Index index; } max_pivot = {0, pivot_index};

    for (Index i = pivot_index; i < n; ++i) {
        const Row &row = mat[i];
        double max_abs_value = max_abs(begin(row) + pivot_index, end(row));  // scaling factor
        if (numeric::isnull(max_abs_value, tol)) return 0;  // null row --> singular matrix
        // update max pivot:
        double scaled_pivot = row[pivot_index]/max_abs_value;
        if (scaled_pivot > max_pivot.value) max_pivot = {scaled_pivot, i};
    }

    std::swap(mat[pivot_index], mat[max_pivot.index]);
    return 1;
}

int lu(double **mat, int n, int perm[], double tol) {
}

