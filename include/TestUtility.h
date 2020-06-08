#ifndef BEAST_INCLUDE_TESTUTILITY_H_
#define BEAST_INCLUDE_TESTUTILITY_H_

#include <cstring>
#include <bitset>
#include "Utility.h"
#include "Trellis.h"
#include "catch.hpp"

void exhaustive_subtrellis_verification(unsigned int n, unsigned int k, int** g, Trellis &trellis, unsigned int w);

#endif //BEAST_INCLUDE_TESTUTILITY_H_
