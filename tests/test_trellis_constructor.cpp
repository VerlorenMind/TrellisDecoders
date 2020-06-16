#include <fstream>
#include <cstdlib>
#include <Utility.h>
#include <Trellis.h>
#include "TestUtility.h"
#include "Simulation.h"

TEST_CASE("CONSTRUCTOR: Can read test matrix from file") {
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);
  int answer[3][6] = {
      {1, 0, 1, 1, 0, 1},
      {1, 0, 1, 0, 1, 0},
      {1, 1, 0, 1, 0, 0},
  };
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      REQUIRE(answer[i][j] == matrix[i][j]);
    }
  }
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form for test matrix is achieved") {
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);
  int answer[3][6] = {
      {1, 1, 0, 1, 0, 0},
      {0, 1, 1, 1, 1, 0},
      {0, 0, 0, 1, 1, 1},
  };
  minspanForm(n, k, matrix);

  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      REQUIRE(answer[i][j] == matrix[i][j]);
    }
  }
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for  test matrix can be found") {
  std::ifstream filename("../tests/test_matrix");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);

  unsigned int *ranges = findRanges(n, k, matrix);

  REQUIRE(ranges[0] == 0);
  REQUIRE(ranges[1] == 5);
  REQUIRE(ranges[2] == 0);
  REQUIRE(ranges[3] == 4);
  REQUIRE(ranges[4] == 0);
  REQUIRE(ranges[5] == 3);

  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
  delete[] ranges;
}

TEST_CASE("CONSTRUCTOR: Can read RM(16,5) matrix from file") {
  std::ifstream filename("../data/reed-muller-16-1");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);

  int answer[5][16] = {
      {1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1},
      {0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1},
      {0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1},
      {0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
  };
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      REQUIRE(answer[i][j] == matrix[i][j]);
    }
  }
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form for RM(16,5) matrix is achieved") {
  std::ifstream filename("../data/reed-muller-16-1");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);

  int answer[5][16] = {
      {1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0},
      {0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0},
      {0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
  };
  minspanForm(n, k, matrix);

  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      REQUIRE(answer[i][j] == matrix[i][j]);
    }
  }
  // WARN("Resulted matrix:\n" << matrixToSstream(k, n, matrix).str());
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for RM(16, 5) matrix can be found") {
  std::ifstream filename("../data/reed-muller-16-1");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);
  minspanForm(n, k, matrix);

  unsigned int *ranges = findRanges(n, k, matrix);

  // WARN("Resulted matrix:\n" << matrixToSstream(k, n, matrix).str());
  REQUIRE(ranges[0] == 0);
  REQUIRE(ranges[1] == 8);
  REQUIRE(ranges[2] == 1);
  REQUIRE(ranges[3] == 14);
  REQUIRE(ranges[4] == 2);
  REQUIRE(ranges[5] == 13);
  REQUIRE(ranges[6] == 3);
  REQUIRE(ranges[7] == 11);
  REQUIRE(ranges[8] == 4);
  REQUIRE(ranges[9] == 15);

  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
  delete[] ranges;
}

TEST_CASE("CONSTRUCTOR: Minspan form for RM(16,5) check matrix is achieved") {
  std::ifstream filename("../data/reed-muller-16-1");
  std::string name;
  unsigned int n, k;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);
  // WARN("Read matrix:\n" << matrixToSstream<int>(k, n, matrix).str());
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
  matrix = readMatrix(filename, n, n - k);
  // WARN("Read matrix:\n" << matrixToSstream<int>(n-k, n, matrix).str());

  int answer[11][16] = {
      {1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1},
  };
  minspanForm(n, n - k, matrix);

  // WARN("Resulted matrix:\n" << matrixToSstream(n-k, n, matrix).str());
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      REQUIRE(answer[i][j] == matrix[i][j]);
    }
  }
  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}

TEST_CASE("TRELLIS: Can reduce a trellis to weight") {
  std::ifstream filename("../tests/test_matrix");
  std::ofstream out;
  std::string name;
  unsigned int n, k;
  unsigned int w = 3;
  std::getline(filename, name);
  filename >> n >> k;
  int **matrix = readMatrix(filename, n, k);
  Trellis trel;
  trel.construct_from_gen_matrix(n, k, matrix);
  out.open("../tests/trellis.gv");
  trel.print_trellis(out);
  out.close();
  system("dot ../tests/trellis.gv -Tpng -o ../tests/before.png");

  trel.reduce_to_weight(w);
  out.open("../tests/trellis.gv");
  trel.print_trellis(out);
  out.close();
  system("dot ../tests/trellis.gv -Tpng -o ../tests/after.png");

  exhaustive_subtrellis_verification(n, k, matrix, trel, w);

  for (unsigned int i = 0; i < k; ++i) {
    delete[] matrix[i];
  }
  delete[] matrix;
}
