#include "SoftDecoder.h"

double SoftDecoder::metric(int x, unsigned int pos) {
  ++op_cmp;
  return (x == alpha[pos] ? 0 : beta[pos]);
}
SoftDecoder::SoftDecoder(unsigned int n, unsigned int k) : n(n), k(k) {
  alpha = new int[n];
  beta = new double[n];
}
SoftDecoder::~SoftDecoder() {
  delete[] alpha;
  delete[] beta;
}
