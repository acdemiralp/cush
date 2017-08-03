#include "catch.hpp"

#include <cush/spherical_harmonics.h>

TEST_CASE("Spherical harmonics coefficient counts are computed.", "[spherical_harmonics]") {
  REQUIRE(cush::coefficient_count(0) == 1);
  REQUIRE(cush::coefficient_count(2) == 9);
  REQUIRE(cush::coefficient_count(4) == 25);
  REQUIRE(cush::coefficient_count(6) == 49);
  REQUIRE(cush::coefficient_count(8) == 81);
}

// Computed by WolframAlpha: SphericalHarmonicY[l, m, theta, phi]
TEST_CASE("Spherical harmonics are computed.", "[spherical_harmonics]") {
  REQUIRE(cush::evaluate(0, 0, M_PI / 2, M_PI / 2) == Approx( 0.2820947918));
  REQUIRE(cush::evaluate(2, 0, M_PI / 2, M_PI / 2) == Approx(-0.3153915652));
  REQUIRE(cush::evaluate(4, 2, M_PI / 2, M_PI / 2) == Approx( 0.3345232718));
  REQUIRE(cush::evaluate(6, 6, M_PI / 2, M_PI / 2) == Approx(-0.4830841135));
  REQUIRE(cush::evaluate(6, 5, M_PI / 2, M_PI / 2) == Approx( 0           ));
  REQUIRE(cush::evaluate(6, 4, M_PI / 2, M_PI / 2) == Approx(-0.3567812628));
  REQUIRE(cush::evaluate(8, 7, M_PI / 2, M_PI / 2) == Approx( 0           ));
  REQUIRE(cush::evaluate(8, 6, M_PI / 2, M_PI / 2) == Approx( 0.3764161087));
  REQUIRE(cush::evaluate(8, 5, M_PI / 2, M_PI / 2) == Approx( 0           ));
}