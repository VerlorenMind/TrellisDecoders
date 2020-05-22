#ifndef BEAST_UTIL_H
#define BEAST_UTIL_H

#include <sstream>
#include <fstream>

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

unsigned int *findRanges(unsigned int n, unsigned int k, int **a);

#endif //BEAST_UTIL_H
