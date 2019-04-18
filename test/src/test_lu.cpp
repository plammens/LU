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
            CHECK(LUDecomposition(test2).decompMatrix() == test2);

            test2 = {{1, 0.5},
                     {1,   1}};
            CHECK(LUDecomposition(test2).decompMatrix() == Matrix{{1, 0.5},
                                                                     {1, 0.5}});

            test2 = {{ 2, 5},
                     {-1, 1}};
            CHECK(LUDecomposition(test2).decompMatrix() == Matrix{{-1, 1},
                                                                     {-2, 7}});

            test2 = {{ 2,   3},
                     {-1, 0.5}};
            CHECK(LUDecomposition(test2).decompMatrix() == Matrix{{-1, 0.5},
                                                                     {-2,   4}});
        }

        SUBCASE("3 by 3") {
            Matrix test3 = {{1, 2, 1},
                            {0, 1, 0},
                            {1, 0, 0}};
            LUDecomposition lu(test3);
            CHECK(lu.decompMatrix() == Matrix{{1., 0., 0.},
                                                 {0., 1., 0.},
                                                 {1., 2., 1.}});
            CHECK((lu.perm().vector() == std::valarray<index_t>{2, 1, 0}).min());
        }

        SUBCASE("singular") {
            CHECK_THROWS_AS(LUDecomposition(Matrix{{0, 0},
                                                   {0, 0}}), SingularMatrixError);
            CHECK_THROWS_AS(LUDecomposition(Matrix{{1, 1},
                                                   {1, 1}}), SingularMatrixError);
            CHECK_THROWS_AS(LUDecomposition(Matrix{{9e-13, 0},
                                                   {0,     1}}), SingularMatrixError);

            Matrix singular = {{1, 0, 0},
                               {0, 1, 1},
                               {1, 1, 1}};

            CHECK_THROWS_AS(LUDecomposition lu(singular);, SingularMatrixError);
        }


    }

    TEST_CASE("array") {
        auto a = newmat(2);
        int perm[2];

        CHECK(not lu(a, 2, perm, numcomp::DEFAULT_TOL));

        a = newmat(2);
        a[0][0] = 1; a[1][1] = 1;
        CHECK(lu(a, 2, perm, numcomp::DEFAULT_TOL) == 1);
        CHECK((a[0][0] == 1 and a[0][1] == 0 and
               a[1][0] == 0 and a[1][1] == 1));
        CHECK((perm[0] == 0 and perm[1] == 1));

        a[0][0] =  2; a[0][1] = 5;
        a[1][0] = -1; a[1][1] = 1;
        CHECK(lu(a, 2, perm, numcomp::DEFAULT_TOL) == -1);
        CHECK((a[0][0] == -1 and a[0][1] == 1 and
               a[1][0] == -2 and a[1][1] == 7));
        CHECK((perm[0] == 1 and perm[1] == 0));
    }

}

