#include "catch.hpp"

#include <cush/wigner.h>

TEST_CASE("Wigner 3J coefficients are computed.", "[wigner]") {
  // Remember parameter units are half-integer.
  REQUIRE(cush::wigner_3j<float>(12, 8, 4, 0, 0, 0) == Approx(0.186989f));
  REQUIRE(cush::wigner_3j<float>( 6, 4, 6, 0, 0, 0) == Approx(0.19518f ));
}
