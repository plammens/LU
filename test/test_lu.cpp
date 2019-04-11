#include <doctest.h>
#include "../src/lu.cpp"


TEST_SUITE("lu") {

    TEST_CASE("max_abs") {
        Row arr = {-1, 2, -3, 1.5};
        CHECK(max_abs(begin(arr), end(arr)) == 3);
    }

}
