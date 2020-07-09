#include <fstream>
#include <ViterbiDecoder.h>
#include "Utility.h"
#include "Simulation.h"
#include "TestUtility.h"
#include "catch.hpp"
const double EPSILON = 0.0000001;

TEST_CASE("VITERBI: Can decode zero word", "[viterbi]") {
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  std::getline(filename, name);
  unsigned int n, k;
  filename >> n >> k;
  int **g;
  g = readMatrix(filename, n, k);
  ViterbiDecoder dec(n, k, g);
  double *x = new double[n];
  int *u = new int[n];
  x[0] = -1.5;
  x[1] = -0.5;
  x[2] = -1.7;
  x[3] = -1.9;
  x[4] = -2.0;
  x[5] = -1.0;

  dec.decode(x, u);

  for (unsigned int i = 0; i < n; ++i) {
    REQUIRE(u[i] == 0);
  }

  delete[] x;
  delete[] u;
  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
}

TEST_CASE("VITERBI: Can decode non-zero word", "[viterbi]") {
  // notice - with this data it should fail
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  std::getline(filename, name);
  int **g, **h;
  unsigned int n, k;
  filename >> n >> k;
  g = readMatrix(filename, n, k);
  h = readMatrix(filename, n, n-k);
  ViterbiDecoder dec(n, k, h);
  double *x = new double[n];
  int *y = new int[n];
  int *u = new int[n];
  u[0] = 0;
  u[1] = 1;
  u[2] = 1;
  u[3] = 0;
  u[4] = 0;
  u[5] = 1;
  x[0] = 0.503536;
  x[1] = 0.423184;
  x[2] = 1.87715;
  x[3] = -0.107348;
  x[4] = -1.7087;
  x[5] = -1.10187;
  // True weight: 1.6053

  std::string codeword = "000000";

  std::ofstream out;
  out.open("../tests/trellis.gv");
  dec.trellis.print_trellis(out);
  out.close();
  system("dot ../tests/trellis.gv -Tpng -o ../tests/trellis.png");

  dec.decode(x, y);

  for (unsigned int j = 0; j < n; ++j) {
    codeword[j] = y[j] ? '1' : '0';
  }
  INFO("Decoded word: " << codeword);
  for (unsigned int i = 0; i < n; ++i) {
    CHECK(y[i] == u[i]);
  }

  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
  delete[] x;
  delete[] u;
  delete[] y;
}


TEST_CASE("VITERBI: Can decode series of random words with minimal noise", "[viterbi]") {
  test_viterbi_decoder(1000, 100, "../tests/test_matrix");
}

TEST_CASE("VITERBI: Can decode series of random words with some noise", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../tests/test_matrix");
}

TEST_CASE("VITERBI: Can decode BCH(7, 4, 3)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-7-4");
}

TEST_CASE("VITERBI: Can decode BCH(7, 1, 5)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-7-1");
}

TEST_CASE("VITERBI: Can decode EBCH(8, 4, 4)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-8-4");
}

TEST_CASE("VITERBI: Can decode BCH(31, 16, 7)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-31-16-7");
}

TEST_CASE("VITERBI: Can decode BCH(31, 21, 5)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-31-21-5");
}

TEST_CASE("VITERBI: Can decode BCH(63, 7, 57)", "[viterbi]") {
  test_viterbi_decoder(1000, 5, "../data/bch-63-7-57");
}

TEST_CASE("VITERBI: Can decode BCH(63, 16, 23)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-63-16-23");
}

/*TEST_CASE("VITERBI: Can decode BCH(63, 30, 13)", "[viterbi]")
{
    test_viterbi_decoder(1000, 1, "../data/bch-63-30-13");
}

TEST_CASE("VITERBI: Can decode BCH(63, 39, 9)", "[viterbi]")
{
    test_viterbi_decoder(1000, 1, "../data/bch-63-39-9");
}*/

TEST_CASE("VITERBI: Can decode RM(16, 5)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/reed-muller-16-1");
}

TEST_CASE("VITERBI: Can decode RM(32, 6)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/reed-muller-32-1");
}

TEST_CASE("VITERBI: Can decode EBCH(32, 16)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-32-16");
}

TEST_CASE("VITERBI: Can decode EBCH(32, 21)", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-32-21");
}

TEST_CASE("VITERBI: Can decode EBCH(32, 16) in SBO", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-32-16-bit-order");
}

TEST_CASE("VITERBI: Can decode EBCH(32, 21) in SBO", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-32-21-bit-order");
}

TEST_CASE("VITERBI: Can decode EBCH(64, 57) in SBO", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-64-57-bit-order");
}

TEST_CASE("VITERBI: Can decode EBCH(64, 51) in SBO", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-64-51-bit-order");
}

TEST_CASE("VITERBI: Can decode EBCH(64, 45) in SBO", "[viterbi]") {
  test_viterbi_decoder(1000, 1, "../data/bch-64-45-bit-order");
}

TEST_CASE("VITERBI: Can decode in reduced trellis", "[viterbi]") {
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  std::getline(filename, name);
  int **g, **h;
  unsigned int n, k;
  filename >> n >> k;
  g = readMatrix(filename, n, k);
  h = readMatrix(filename, n, n-k);
  ViterbiDecoder dec(n, k, h);
  double *x = new double[n];
  int *y = new int[n];
  int *u = new int[n];
  int *u_around = new int[n];
  u[0] = 0;
  u[1] = 1;
  u[2] = 1;
  u[3] = 1;
  u[4] = 1;
  u[5] = 0;
  u_around[0] = 0;
  u_around[1] = 0;
  u_around[2] = 0;
  u_around[3] = 1;
  u_around[4] = 1;
  u_around[5] = 1;
  x[0] = -0.503536;
  x[1] = 0.423184;
  x[2] = 1.87715;
  x[3] = -0.107348;
  x[4] = 1.7087;
  x[5] = -1.10187;
  // True weight: -0.107348

  std::string codeword = "000000";

  dec.trellis.reduce_to_weight(3);
  std::ofstream out;
  out.open("../tests/trellis.gv");
  dec.trellis.print_trellis(out);
  out.close();
  system("dot ../tests/trellis.gv -Tpng -o ../tests/trellis.png");

  dec.decode_around(x, y, u_around);

  for (unsigned int j = 0; j < n; ++j) {
    codeword[j] = y[j] ? '1' : '0';
  }
  INFO("Decoded word: " << codeword);
  for (unsigned int i = 0; i < n; ++i) {
    CHECK(y[i] == u[i]);
  }

  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
  delete[] x;
  delete[] u;
  delete[] y;
}
