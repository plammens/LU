//
// Created by Paolo on 16/05/2019.
//

#include <doctest.h>
#include <iostream>

#include "debug.h"
#include "inverse.h"
#include "norm.h"
#include "extra.h"


TEST_SUITE ("inverse") {

    TEST_CASE ("simple") {
        Matrix mat = {{1, 2},
                      {3, 4}};
        Matrix inv = inverse(mat);
        CHECK(normInf(mat*inv - identity(2)) < 1e-12);
    }

    TEST_CASE ("rcond") {
        Matrix Id = identity(10);
        CHECK(conditionNumber(Id, NormType::L1) == doctest::Approx(1.0));
        CHECK(conditionNumber(Id, NormType::Inf) == doctest::Approx(1.0));
    }

    TEST_CASE ("determinant") {
        Matrix A(10);
        CHECK(determinant(Matrix(10)) == 0.0);
        CHECK(determinant(Matrix{{9e-13, 0},
                                 {0,     1}}) == 0.0);
        CHECK(determinant(Matrix{{1, -1, 3},
                                 {4, 5, 2},
                                 {3, 1, 2}}) == doctest::Approx(-23.0));
    }

}
