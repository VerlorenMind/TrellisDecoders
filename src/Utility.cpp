#include <iostream>
#include <Exceptions.h>
#include "Utility.h"

int **readMatrix(std::ifstream &input, unsigned int n, unsigned int k) {
  // Check n and k before allocating ungodly amounts of memory due to corrupted values
  unsigned long long trying_to_alloc = k*n*sizeof(int);
  if(trying_to_alloc > (uint32_t(1) << 30)) {
    throw InvalidFileException();
  }
  int **matrix = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    matrix[i] = new int[n];
  }
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      int temp;
      if(!input.good()) {
        throw InvalidFileException();
      }
      input >> temp;
      if(temp != 0 && temp != 1) {
        throw InvalidFileException();
      }
      matrix[i][j] = temp;
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
    for (unsigned int j = fixed_rows; j < k; ++j) {
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

void minspanForm(unsigned int n, unsigned int k, int **a, std::vector<unsigned> &row_start, std::vector<unsigned> &row_end) {
  std::fill(row_start.begin(), row_start.end(), ~0);
  std::fill(row_end.begin(), row_end.end(), ~0);
  unsigned int num; // number of rows with 1s in i'th column
  unsigned int fixed_rows = 0; // number of rows on top of the matrix that were transformed
  unsigned int *rows = new unsigned int[k]; // rows with 1s in i'th column
  // Left side
  int i = 0;
  while (fixed_rows < k) {
    num = 0;
    for (unsigned int j = fixed_rows; j < k; ++j) {
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
      row_start[i-1] = fixed_rows;
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
      row_end[i+1] = rows[0];
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

unsigned int op_cmp;

unsigned int partition(int lo, int hi, double *x, std::vector<int> &permutation) {
  double piv = x[hi];
  unsigned int i = lo;
  for (unsigned int j = lo; j <= hi; ++j) {
    if (x[j] > piv) {
      std::swap(x[i], x[j]);
      std::swap(permutation[i], permutation[j]);
      ++i;
    }
    ++op_cmp;
  }
  std::swap(x[i], x[hi]);
  std::swap(permutation[i], permutation[hi]);
  return i;
}

void quicksort(int lo, int hi, double *x, std::vector<int> &permutation) {
  if (lo < hi) {
    unsigned int p = partition(lo, hi, x, permutation);
    quicksort(lo, p - 1, x, permutation);
    quicksort(p + 1, hi, x, permutation);
  }
}

unsigned int quicksort(unsigned int n, double *x, std::vector<int> &permutation) {
  op_cmp = 0;
  permutation.resize(n);
  for (unsigned int i = 0; i < n; ++i) {
    permutation[i] = i;
  }
  quicksort(0, n - 1, x, permutation);
  return op_cmp;
}

void echelonForm(unsigned int n, unsigned int k, int **x) {
  unsigned int r = 0;
  unsigned int c = 0;
  while (r < k && c < n) {
    int index = -1;
    for (unsigned int i = r; i < k; ++i) {
      if (x[i][c]) {
        index = i;
        break;
      }
    }
    if (index == -1) {
      ++c;
    } else {
      if (index != r) {
        for (unsigned int i = 0; i < n; ++i) {
          std::swap(x[index][i], x[r][i]);
        }
      }
      for (unsigned int i = r + 1; i < k; ++i) {
        if (x[i][c]) {
          for (unsigned int j = c; j < n; ++j) {
            x[i][j] ^= x[r][j];
          }
        }
      }
      ++r;
      ++c;
    }
  }
}

void comb_search(unsigned int n,
                 unsigned int k,
                 int offset,
                 std::vector<int> &comb,
                 std::vector<std::vector<int>> &combs) {
  if (k == 0) {
    combs.push_back(comb);
    return;
  } else {
    for (unsigned int i = offset; i < n; ++i) {
      comb.push_back(i);
      comb_search(n, k - 1, i + 1, comb, combs);
      comb.pop_back();
    }
  }
}

std::vector<std::vector<int>> combinations(unsigned int n, unsigned int k) {
  std::vector<std::vector<int>> combs(0);
  std::vector<int> comb(0);
  comb_search(n, k, 0, comb, combs);
  return combs;
}

std::vector<int> weight_profile(int n, int k, int **g) {
  std::vector<int> weight_profile;
  weight_profile.resize(0);
  int word_weight = 0;
  int *codeword = new int[n];
  for (uint64_t infoword = 0; infoword < (uint64_t(1) << k); ++infoword) {
    if (infoword == 0) {
      memset(codeword, 0, n * sizeof(int));
    } else {
      uint64_t bit_diff = infoword ^(infoword - 1);
      for (unsigned int l = 0; l < k; ++l) {
        if (bit_diff & (uint64_t(1) << l)) {
          for (unsigned int j = 0; j < n; ++j) {
            codeword[j] ^= g[l][j];
          }
        }
      }
    }
    word_weight = 0;
    for (unsigned int i = 0; i < n; ++i) {
      word_weight += codeword[i];
    }
    auto iter = weight_profile.begin();
    while (iter != weight_profile.end() && *iter < word_weight) {
      ++iter;
    }
    if (iter == weight_profile.end()) {
      weight_profile.push_back(word_weight);
    } else if (*iter != word_weight) {
      weight_profile.insert(iter, word_weight);
    }
  }
  delete[] codeword;
  return weight_profile;
}
