#include "catch.hpp"

#include <cush/clebsch_gordan.h>

TEST_CASE("Clebsch Gordan coefficients are computed.", "[clebsch_gordan]") {
  REQUIRE(cush::clebsch_gordan<float>(5, 4, 1, 0, 0, 0) == Approx( 0.3892494f));
  REQUIRE(cush::clebsch_gordan<float>(3, 2, 3, 0, 0, 0) == Approx(-0.5163977f));
}
