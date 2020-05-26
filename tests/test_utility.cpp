#include <cstring>
#include "catch.hpp"
#include "Utility.h"

TEST_CASE("UTIL: Test quicksort") {
  std::vector<int> perm = {8, 4, 7, 3, 6, 2, 5, 1, 0};
  double x[9] = {-100, -10, 0, 10, 100, -5, 5, 50, 1000};
  std::vector<int> perm1;

  unsigned int ops = quicksort(9, x, perm1);

  REQUIRE(ops > 0);
  for(unsigned int i=0; i<9; ++i) {
    REQUIRE(perm[i] == perm1[i]);
    if(i > 0) {
      REQUIRE(x[i] < x[i-1]);
    }
  }
}
