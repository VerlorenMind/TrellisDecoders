#ifndef BEAST_INCLUDE_VITERBIDECODER_H_
#define BEAST_INCLUDE_VITERBIDECODER_H_

#include "TrellisDecoder.h"
#include "Trellis.h"
class ViterbiDecoder : public SoftDecoder {
 private:
  double *prev_layer, *cur_layer;
  uint64_t *prev_layer_path, *cur_layer_path;
  double decode(int *u);
 public:
  Trellis trellis;
  ViterbiDecoder(unsigned int n, unsigned int k, int **g);
  ~ViterbiDecoder();
  double decode(const double *x, int *u) override;
  double decode_around(const double *x, int *u, const int *c);
};
#endif //BEAST_INCLUDE_VITERBIDECODER_H_
