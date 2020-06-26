#ifndef BEAST_INCLUDE_TESTUTILITY_H_
#define BEAST_INCLUDE_TESTUTILITY_H_

#include <cstring>
#include <bitset>
#include "Utility.h"
#include "Trellis.h"
#include "Simulation.h"
#include <OrderedStatisticsDecoder.h>
#include "catch.hpp"

void exhaustive_subtrellis_verification(unsigned int n, unsigned int k, int** g, Trellis &trellis, unsigned int w, int **h=nullptr);

void test_beast_decoder(int tests, double stn, double delta, const std::string &filename);
void test_bsd_decoder(int tests, double stn, std::string filename);
void test_KTKL_decoder(int tests, double stn, int w, int buf_size, const std::string &filename);
void test_osd_decoder(int tests, double stn, int w, std::string filename);
void test_viterbi_decoder(int tests, double stn, const std::string &filename);
#endif //BEAST_INCLUDE_TESTUTILITY_H_
