#include <iostream>
#include "Utility.h"

int **readMatrix(std::ifstream &input, unsigned int n, unsigned int k) {
  int **matrix = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    matrix[i] = new int[n];
  }
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      input >> matrix[i][j];
    }
  }
  return matrix;
}

void minspanForm(unsigned int n, unsigned int k, int **a) {
  unsigned int num; // number of rows with 1s in i'th column
  unsigned int fixed_rows = 0; // number of rows on top of the matrix that were transformed
  unsigned int *rows = new unsigned int[k]; // rows with 1s in i'th column
  std::string tempmatr;
  // Left side
  int i = 0;
  while (fixed_rows < k) {
    tempmatr = matrixToSstream(k, n, a).str();
    num = 0;
    for (int j = fixed_rows; j < k; ++j) {
      if (a[j][i]) {
        rows[num++] = (unsigned int) j;
      }
    }
    ++i;
    if (num == 0) {
      continue;
    } else {
      // If the first row after fixed rows does not contain 1, swap it with the first row that does
      if (rows[0] != fixed_rows) {
        for (unsigned int l = 0; l < n; ++l) {
          std::swap(a[rows[0]][l], a[fixed_rows][l]);
        }
        rows[0] = fixed_rows;
      }

      ++fixed_rows;

      // Adding fixed row to the unfixed rest of the matrix
      for (unsigned int l = 1; l < num; ++l) {
        for (unsigned int j = 0; j < n; ++j) {
          a[rows[l]][j] ^= a[rows[0]][j];
        }
      }
    }
  }
  // Right side
  // Same stuff as above, but with different indices and without swapping rows
  fixed_rows = 0;
  int *fixed_nums = new int[k];
  i = n - 1;
  while (fixed_rows < k) {
    tempmatr = matrixToSstream(k, n, a).str();
    num = 0;
    for (int j = k - 1; j >= 0; --j) {
      bool flag = true;
      for (unsigned int l = 0; l < fixed_rows; ++l) {
        if (fixed_nums[l] == j) {
          flag = false;
          break;
        }
      }
      if (flag) {
        if (a[j][i]) {
          rows[num++] = (unsigned int) j;
        }
      }
    }
    --i;
    if (num == 0) {
      continue;
    } else {
      fixed_nums[fixed_rows++] = rows[0];
      for (unsigned int l = 1; l < num; ++l) {
        for (unsigned int j = 0; j < n; ++j) {
          a[rows[l]][j] ^= a[rows[0]][j];
        }
      }
    }
  }
  delete[] rows;
  delete[] fixed_nums;
}

unsigned int *findRanges(unsigned int n, unsigned int k, int **a) {
  // Finding active ranges for each row (range between the first non-zero element and the last one)
  unsigned int *ranges = new unsigned int[2 * k];
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      if (a[i][j]) {
        ranges[2 * i] = j;
        break;
      }
    }
    for (int j = n - 1; j >= 0; --j) {
      if (a[i][j]) {
        ranges[2 * i + 1] = (unsigned int) j;
        break;
      }
    }
  }
  return ranges;
}


