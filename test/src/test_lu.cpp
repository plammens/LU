#include <doctest.h>
#include <iostream>
#include <numeric>
#include "debug.h"


//TEST_SUITE ("lu") {
//
//    TEST_CASE ("max_abs") {
//        Row_ arr = {-1, 2, -3, 1.5};
//        CHECK(max_abs(begin(arr), end(arr)) == 3);
//    }
//
//    TEST_CASE ("scaled partial pivoting") {
//        Matrix_ mat = {{1, 2},
//                      {3, 4}};
//        REQUIRE(scaledPartialPivoting(mat, 0));
//        CHECK(matEqual(mat, Matrix_{{3, 4},
//                                   {1, 2}}));
//
//        mat = {{10, 100},
//               {5,  5}};
//        CHECK(scaledPartialPivoting(mat, 0));
//        CHECK(matEqual(mat, Matrix_{{5,  5},
//                                   {10, 100}}));
//
//        SUBCASE("singular") {
//            mat = {{0, 0},
//                   {0, 0}};
//            REQUIRE(not scaledPartialPivoting(mat, 0));
//
//            mat = {{1,       1},
//                   {0.9e-12, 0.9e-12}};
//            CHECK(not scaledPartialPivoting(mat, 0));
//        }
//    }
//
//}


TEST_SUITE ("matrix") {

    TEST_CASE ("ctor") {
        Matrix mat(3);
        REQUIRE(mat.size() == 3);
        REQUIRE(mat.data().size() == 9);
        CHECK(mat == Matrix({{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0}}));

        SUBCASE("init_list") {
            CHECK_THROWS_AS(Matrix({{1, 2}, {1}}), std::invalid_argument);
            Matrix mat = {{1, 2}, {3, 4}};
            REQUIRE(mat.size() == 2);
            REQUIRE(mat == mat);
        }
    }

    TEST_CASE("subscript") {
        Matrix mat = {{1, 2},
                      {3, 4}};

        CHECK_THROWS_AS(mat(0, 2), std::out_of_range);
        CHECK_THROWS_AS(mat(2, 0), std::out_of_range);
        CHECK_THROWS_AS(mat(2, 2), std::out_of_range);
        CHECK_THROWS_AS(mat(-1, 0), std::out_of_range);

        bool correctElements = mat(0, 0) == 1 and mat(0, 1) == 2 and
                               mat(1, 0) == 3 and mat(1, 1) == 4;
        CHECK(correctElements);
    }

    TEST_CASE("iterators") {
        Matrix mat = {{1,    2},
                      {3.14, 4}};

        double sum = 0;
        for (auto it = mat.begin(); it < mat.end(); ++it)
            sum += *it;

        double expectedSum = mat.data().sum();
        CHECK(sum == expectedSum);
        CHECK(std::accumulate(mat.begin(), mat.end(), 0.0) == expectedSum);
    }

}


TEST_SUITE("triangular") {
    Matrix null3 = {{0, 0, 0},
                    {0, 0, 0},
                    {0, 0, 0}};

    TEST_CASE("constructor") {

        SUBCASE("default") {
            LowerTriangular lower(3);
            UpperTriangular upper(3);
            REQUIRE(lower.size() == 3);
            CHECK(lower == null3);
            REQUIRE(upper.size() == 3);
            CHECK(upper == null3);
            CHECK(lower == upper);
        }

        SUBCASE("lower") {
            CHECK_THROWS_AS(LowerTriangular({{1, 2}, {3, 4}}), std::invalid_argument);
            REQUIRE_NOTHROW(LowerTriangular({{1}, {2, 3}}));
            REQUIRE_NOTHROW(LowerTriangular({{1, 0}, {2, 3}}));

            const LowerTriangular mat{{1},
                                      {2, 3.14}};

            CHECK(mat == Matrix{{1, 0},
                                {2, 3.14}});
            CHECK(mat(0, 1) == 0);
        }

        SUBCASE("upper") {
            CHECK_THROWS_AS(UpperTriangular({{1, 2}, {3, 4}}), std::invalid_argument);
            REQUIRE_NOTHROW(UpperTriangular({{1, 2}, {0, 3}}));

            const UpperTriangular mat{{1, 2},
                                      {0, 3.14}};

            CHECK(mat == Matrix{{1, 2},
                                {0, 3.14}});

            CHECK(mat(1, 0) == 0);
        }

    }

}
