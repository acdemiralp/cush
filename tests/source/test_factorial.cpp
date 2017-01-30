#include <catch.hpp>

#include <factorial.h>

TEST_CASE("2- Factorials are computed.", "[factorial]") {
  REQUIRE(cush::factorial<float>(0)  == 1);
  REQUIRE(cush::factorial<float>(1)  == 1);
  REQUIRE(cush::factorial<float>(2)  == 2);
  REQUIRE(cush::factorial<float>(3)  == 6);
  REQUIRE(cush::factorial<float>(10) == 3628800);
}
