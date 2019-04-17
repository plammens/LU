#include <doctest.h>
#include <iostream>
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
        std::cout << mat << std::endl;
        CHECK(correctElements);
    }

    TEST_CASE("iterators") {
        Matrix mat = {{1,    2},
                      {3.14, 4}};

        double sum = 0;
        for (auto it = mat.begin(); it < mat.end(); ++it)
            sum += *it;

        double expectedSum = mat.data().sum().sum();
        CHECK(sum == expectedSum);
        CHECK(std::accumulate(mat.begin(), mat.end(), 0.0) == expectedSum);
    }

}


TEST_SUITE("lu") {

    TEST_CASE("simple") {
        Matrix test = {{1, 0},
                       {0, 1}};
        CHECK(LUDecomposition(test).getDecompMatrix() == test);

        test = {{1, 0.5},
                {1,   1}};
        CHECK(LUDecomposition(test).getDecompMatrix() == Matrix{{1, 0.5},
                                                                {1, 0.5}});

        Matrix test3 = {{1, 2, 1},
                        {0, 1, 0},
                        {1, 0, 0}};
        std::cout << LUDecomposition(test3).getDecompMatrix() << std::endl;
    }

}

