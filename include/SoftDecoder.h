#ifndef BEAST_ISOFTDECODER_H
#define BEAST_ISOFTDECODER_H

#include "DecoderID.h"

class SoftDecoder {
 protected:
  DecoderID id = DecoderID::ERROR;
  unsigned int n;
  unsigned int k;
  int *alpha;
  double *beta;
 public:
  long op_cmp, op_add;
  double metric(int x, unsigned int pos);
  SoftDecoder(unsigned int n, unsigned int k);
  ~SoftDecoder();
  virtual double decode(double *x, int *u) = 0;
  DecoderID get_id() { return id; };
};
#endif //BEAST_ISOFTDECODER_H
