#include <doctest.h>
#include "debug.h"

#include "sistema.h"
#include "oldies.h"


TEST_SUITE("sistema") {

    TEST_CASE("simple") {

        SUBCASE("identity") {
            Matrix id = {{1, 0},
                         {0, 1}};

            SolveResult result = solve(id, {3.14, 2.78});
            CHECK(result);
            CHECK((result.solution() == Vector{3.14, 2.78}).min());
        }

        SUBCASE("singular") {
            Matrix mat(3);

            SolveResult result = solve(mat, {});
            CHECK(not result);
            CHECK_THROWS_AS(result.solution(), std::invalid_argument);
        }

        SUBCASE("numerical") {
            Matrix mat = {{1,  1e10},
                          {1, 1e-10}};
            auto result = solve(mat, {1e10, 1});
            CHECK(result);
            CHECK(result.solution()[0] == doctest::Approx(1e10/(1e10 + 1)));
            CHECK(result.solution()[1] == doctest::Approx(1e10/(1e10 + 1)));

            mat = {{-1e-20, 1},
                   {    2, 1}};
            result = solve(mat, {1, 0});
            CHECK(result);
            CHECK(result.solution()[0] == doctest::Approx(-0.5).epsilon(1e-20));
            CHECK(result.solution()[1] == doctest::Approx(1).epsilon(1e-20));
        }

    }

    TEST_CASE("C-style") {

        auto a = newmat(2);
        double b[2] = {1e10, 1};
        double x[2];

        CHECK(not sistema(a, x, b, 2, numcomp::DEFAULT_TOL));

        a = newmat(2);
        a[0][0] = 1; a[0][1] = 1e10;
        a[1][0] = 1; a[1][1] = 1e-10;
        CHECK(sistema(a, x, b, 2, numcomp::DEFAULT_TOL) == -1);
        CHECK(x[0] == doctest::Approx(1e10/(1e10 + 1)));
        CHECK(x[1] == doctest::Approx(1e10/(1e10 + 1)));

        freemat(a, 2);

    }

}


