#ifndef BEAST_TESTUTIL_H
#define BEAST_TESTUTIL_H

#include <iostream>
#include <stack>
#include <random>
#include <chrono>
#include "catch.hpp"
#include "BeastDecoder.h"
#include "BSDDecoder.h"

void generate_vector(unsigned int size, unsigned int* vector);

void test_decoder(unsigned int tests, double stn, double delta, const char* filename, std::string decoder);

#endif //BEAST_TESTUTIL_H
