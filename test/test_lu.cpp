#include <doctest.h>
#include "../src/lu.cpp"


bool matEqual(const Matrix &lhs, const Matrix &rhs) {
    for (auto r1 = begin(lhs), r2 = begin(rhs); r1 != end(lhs); ++r1, ++r2)
        for (auto it1 = begin(*r1), it2 = begin(*r2); it1 != end(*r1); ++it1, ++it2)
            if (not(*it1 == *it2)) return false;
    return true;
}


TEST_SUITE ("lu") {

    TEST_CASE ("max_abs") {
        Row arr = {-1, 2, -3, 1.5};
        CHECK(max_abs(begin(arr), end(arr)) == 3);
    }

    TEST_CASE ("scaled partial pivoting") {
        Matrix mat = {{1, 2},
                      {3, 4}};
        REQUIRE(scaled_partial_pivoting(mat, 0));
        CHECK(matEqual(mat, Matrix{{3, 4},
                                   {1, 2}}));

        mat = {{10, 100},
               {5,  5}};
        CHECK(scaled_partial_pivoting(mat, 0));
        CHECK(matEqual(mat, Matrix{{5,  5},
                                   {10, 100}}));

        SUBCASE("singular") {
            mat = {{0, 0},
                   {0, 0}};
            REQUIRE(not scaled_partial_pivoting(mat, 0));

            mat = {{1,       1},
                   {0.9e-12, 0.9e-12}};
            CHECK(not scaled_partial_pivoting(mat, 0));
        }
    }

}
