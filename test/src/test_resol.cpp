#include <doctest.h>
#include "debug.h"

#define RESOL_DECLARATIONS_ONLY
#include "resol.h"
#include "oldies.h"


TEST_SUITE("resol") {

    TEST_CASE("simple") {
        SUBCASE("identity") {
            Matrix id = {{1, 0},
                          {0, 1}};
            LUDecomposition luDecomp(id);

            CHECK((solveLU(luDecomp, {3.14, 2.78}) == Vector{3.14, 2.78}).min());
            CHECK((solveLU(luDecomp, {0, 0}) == Vector{0, 0}).min());
        }

        SUBCASE("1") {
            Matrix mat = {{1, 1},
                          {0, 1}};

            LUDecomposition luDecomp(mat);
            CHECK((solveLU(luDecomp, {2, 3}) == Vector{-1, 3}).min());
        }

        SUBCASE("2") {
            Matrix mat = {{0, 1},
                          {1, 1}};

            LUDecomposition luDecomp(mat);
            CHECK((solveLU(luDecomp, {2, 3}) == Vector{1, 2}).min());
        }

        SUBCASE("3") {
            Matrix mat = {{3,    5},
                          {-2, 1.4}};
            LUDecomposition luDecomp(mat);

            Vector x = solveLU(luDecomp, {1, 0});
            CHECK(x[0] == 7./71);
            CHECK(x[1] == 10./71);
        }
    }

    TEST_CASE("array") {
        double **a = newmat(2);
        a[0][0] = 1; a[0][1] = 1; a[1][1] = 1;
        int perm[2];
        lu(a, 2, perm, numcomp::DEFAULT_TOL);

        double b[] = {2, 3};
        double x[2];
        resol(a, x, b, 2, perm);
        CHECK((x[0] == -1 and x[1] == 3));

        freemat(a, 2);
    }

}

