#include <catch.hpp>

#include <legendre.h>

TEST_CASE("4- Associated Legendre polynomials are computed.", "[legendre]") {
  REQUIRE(cush::associated_legendre(0 , 0, 0.5) == Approx(1             ));
  REQUIRE(cush::associated_legendre(2 , 0, 0.5) == Approx(-0.125        ));
  REQUIRE(cush::associated_legendre(2 , 2, 0.5) == Approx(2.25          ));
  REQUIRE(cush::associated_legendre(6 , 4, 0.5) == Approx(465.1171875   ));
  REQUIRE(cush::associated_legendre(10, 6, 0.5) == Approx(-82397.3785400));
}
