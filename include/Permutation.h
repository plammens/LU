//
// Created by Paolo on 14/05/2019.
//

#ifndef LU_PERMUTATION_H
#define LU_PERMUTATION_H

#include <valarray>
#include "aliases.h"


/// Representation of a permutation
class Permutation {
public:
    typedef std::valarray<index_t> Vector;

    explicit Permutation(size_t n);
    Permutation() = default;

    void permute(index_t a, index_t b);
    inline const Vector &vector() const { return _vec; }
    inline bool parity() const { return _parity; }  ///< 0 if even, 1 if odd

private:
    Vector _vec;
    bool _parity = false;
};


#endif //LU_PERMUTATION_H
