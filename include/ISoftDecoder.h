#ifndef BEAST_ISOFTDECODER_H
#define BEAST_ISOFTDECODER_H

#include "DecoderID.h"

class ISoftDecoder {
 protected:
  DecoderID id = DecoderID::ERROR;
 public:
  long op_cmp, op_add;
  ISoftDecoder() = default;
  virtual ~ISoftDecoder() = default;
  virtual double decode(double *x, int *u) = 0;
  DecoderID get_id() { return id; };
};
#endif //BEAST_ISOFTDECODER_H
