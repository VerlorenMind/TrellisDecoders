#ifndef BEAST_INCLUDE_VITERBIDECODER_H_
#define BEAST_INCLUDE_VITERBIDECODER_H_

#include "TrellisDecoder.h"
#include "Trellis.h"
class ViterbiDecoder : public SoftDecoder {
 private:
  Trellis trellis;
  double *prev_layer, *cur_layer;
  uint64_t *prev_layer_path, *cur_layer_path;
 public:
  ViterbiDecoder(unsigned int n, unsigned int k, int **g);
  virtual double decode(double *x, int *u) override;
};
#endif //BEAST_INCLUDE_VITERBIDECODER_H_
