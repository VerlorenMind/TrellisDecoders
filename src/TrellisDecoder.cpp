#include "TrellisDecoder.h"

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, std::ifstream &filename) : SoftDecoder(n, k) {
  int **checkmatrix = readMatrix(filename, n, k);
  init(n, k, checkmatrix);
  for (unsigned int i = 0; i < k; ++i) {
    delete[] checkmatrix[i];
  }
  delete[] checkmatrix;
}

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, int **checkmatrix) : SoftDecoder(n, k) {
  init(n, k, checkmatrix);
}

TrellisDecoder::~TrellisDecoder() {
  for (unsigned int i = 0; i < k; ++i) {
    delete[] h[i];
  }
  delete[] h;
  delete[] ranges;
  delete[] trellis_profile;
}

void TrellisDecoder::init(unsigned int n, unsigned int k, int **checkmatrix) {
// Transforming matrix to minspan form and finding active ranges of it's rows
  h = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    h[i] = new int[n];
    memcpy(h[i], checkmatrix[i], n * sizeof(int));
  }
  minspanForm(n, k, h);
  ranges = findRanges(n, k, h);
  trellis_profile = new uint64_t[n + 1];
  unsigned pow;
  trellis_profile[0] = 1;
  trellis_size = 1;
  for (unsigned i = 0; i < n; ++i) {
    pow = 0;
    for (unsigned int j = 0; j < k; ++j) {
      if (i >= ranges[2 * j] && i < ranges[2 * j + 1]) {
        ++pow;
      }
    }
    trellis_profile[i + 1] = uint64_t(1) << pow;
    trellis_size += trellis_profile[i + 1];
  }
}
