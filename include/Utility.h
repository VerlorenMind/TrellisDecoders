#ifndef BEAST_UTIL_H
#define BEAST_UTIL_H

#include <sstream>
#include <fstream>
#include <vector>

template<class T>
std::stringstream arrayToSstream(unsigned int size, T *a) {
  std::stringstream result;
  for (unsigned int i = 0; i < size; ++i) {
    result << a[i] << " ";
  }
  return result;
}

template<class T>
std::stringstream matrixToSstream(unsigned int k, unsigned int n, T **a) {
  std::stringstream result;
  for (unsigned int i = 0; i < k; ++i) {
    for (unsigned int j = 0; j < n; ++j) {
      result << a[i][j] << " ";
    }
    result << "\n";
  }
  return result;
}

int **readMatrix(std::ifstream &input, unsigned int n, unsigned int k);

void minspanForm(unsigned int n, unsigned int k, int **a);
void minspanForm(unsigned int n, unsigned int k, int **a, std::vector<unsigned> &row_start, std::vector<unsigned> &row_end);

unsigned int *findRanges(unsigned int n, unsigned int k, int **a);

unsigned int quicksort(unsigned int n, double *x, std::vector<int> &permutation);

void echelonForm(unsigned int n, unsigned int k, int **x);

std::vector<std::vector<int>> combinations(unsigned int n, unsigned int k);

#endif //BEAST_UTIL_H
