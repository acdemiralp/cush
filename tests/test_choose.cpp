#include "catch.hpp"

#include <cush/choose.h>

TEST_CASE("Binomial coefficients are computed.", "[choose]") {
  REQUIRE(cush::choose( 0,  0) == 1);
  REQUIRE(cush::choose( 0,  2) == 0);
  REQUIRE(cush::choose( 4,  2) == 6);
  REQUIRE(cush::choose( 7,  5) == 21);
  REQUIRE(cush::choose(24, 16) == 735471);
}
