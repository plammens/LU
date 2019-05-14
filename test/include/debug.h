#ifndef LU_DEBUG_H
#define LU_DEBUG_H

#include <iosfwd>

#include "Matrix.h"

inline
bool operator==(const Matrix &lhs, const Matrix &rhs) {
    return lhs.size() == rhs.size()
           and std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

inline
std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
    std::size_t n = mat.size();
    auto it = mat.begin();

    os << '{';
    for (unsigned i = 0; i < n; ++i) {
        if (i) os << ' ';
        os << '{';
        for (unsigned j = 0; j < n; ++j) {
            if (j) os << ' ';
            os << *it;
            ++it;
        }
        os << '}';
    }
    os << '}';

    return os;
}

#endif //LU_DEBUG_H
