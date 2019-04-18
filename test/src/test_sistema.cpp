#include <doctest.h>
#include "debug.h"

#define SISTEMA_DECLARATIONS_ONLY
#include "sistema.cpp"


TEST_SUITE("sistema") {

    TEST_CASE("simple") {

        SUBCASE("identity") {
            Matrix id = {{1, 0},
                         {0, 1}};

            SolveResult result = solve(id, {3.14, 2.78});
            CHECK(result.success);
            CHECK((result.solution == Vector{3.14, 2.78}).min());
        }

        SUBCASE("singular") {
            Matrix mat(3);

            SolveResult result = solve(mat, {});
            CHECK(not result.success);
            CHECK(result.solution.size() == 0);
        }

        SUBCASE("numerical") {
            Matrix mat = {{1,  1e10},
                          {1, 1e-10}};
            auto result = solve(mat, {1e10, 1});
            CHECK(result.success);
            CHECK(result.solution[0] == doctest::Approx(1e10/(1e10 + 1)));
            CHECK(result.solution[1] == doctest::Approx(1e10/(1e10 + 1)));

            mat = {{-1e-20, 1},
                   {    2, 1}};
            result = solve(mat, {1, 0});
            CHECK(result.success);
            CHECK(result.solution[0] == doctest::Approx(-0.5).epsilon(1e-20));
            CHECK(result.solution[1] == doctest::Approx(1).epsilon(1e-20));
        }

    }

}


