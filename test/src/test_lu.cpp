#include <doctest.h>
#include <numeric>
#include "debug.h"


TEST_SUITE ("matrix") {

    TEST_CASE ("ctor") {
        Matrix mat(3);
        REQUIRE(mat.size() == 3);
        REQUIRE(mat.data().size() == 3);
        CHECK(mat == Matrix({{0, 0, 0},
                             {0, 0, 0},
                             {0, 0, 0}}));

        SUBCASE("init_list") {
            CHECK_THROWS_AS(Matrix({{1, 2}, {1}}), std::invalid_argument);
            Matrix mat2 = {{1, 2}, {3, 4}};
            CHECK(mat2.size() == 2);
            CHECK(mat2 == mat2);
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

        double sum = 0.;
        for (auto it = mat.begin(); it < mat.end(); ++it)
            sum += *it;

        double expectedSum = mat.data().sum().sum();
        CHECK(sum == expectedSum);
        CHECK(std::accumulate(mat.begin(), mat.end(), 0.0) == expectedSum);
    }

}


TEST_SUITE("lu") {

    TEST_CASE("simple") {
        SUBCASE("2 by 2") {
            Matrix test2 = {{1, 0},
                           {0, 1}};
            CHECK(LUDecomposition(test2).getDecompMatrix() == test2);

            test2 = {{1, 0.5},
                     {1,   1}};
            CHECK(LUDecomposition(test2).getDecompMatrix() == Matrix{{1, 0.5},
                                                                     {1, 0.5}});

            test2 = {{ 2, 5},
                     {-1, 1}};
            CHECK(LUDecomposition(test2).getDecompMatrix() == Matrix{{-1, 1},
                                                                     {-2, 7}});

            test2 = {{ 2,   3},
                     {-1, 0.5}};
            CHECK(LUDecomposition(test2).getDecompMatrix() == Matrix{{-1, 0.5},
                                                                     {-2,   4}});
        }

        SUBCASE("singular") {
            CHECK_THROWS_AS(LUDecomposition(Matrix{{0, 0}, {0, 0}}), SingularMatrixError);
            CHECK_THROWS_AS(LUDecomposition(Matrix{{1, 1}, {1, 1}}), SingularMatrixError);
        }

        SUBCASE("3 by 3") {
            Matrix test3 = {{1, 2, 1},
                            {0, 1, 0},
                            {1, 0, 0}};
            CHECK(LUDecomposition(test3).getDecompMatrix() == Matrix{{1., 0., 0.},
                                                                     {0., 1., 0.},
                                                                     {1., 2., 1.}});

            test3 = {{0, 0, 0},
                     {0, 0, 0},
                     {0, 0, 0}};
        }
    }

}

