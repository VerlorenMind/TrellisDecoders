#ifndef BEAST_TRELLISDECODER_H
#define BEAST_TRELLISDECODER_H

#include <fstream>
#include <cstring>

#include "Utility.h"
#include "SoftDecoder.h"

enum Tree {
  NIL,
  FWD,
  BKW,
};

struct Node {
  Tree tree;
  unsigned layer;
  uint64_t number;
  mutable double metric;
  mutable uint64_t path;
  mutable bool path_avalaible[2];
};

class TrellisDecoder : public SoftDecoder {
 protected:
  int **h;
  unsigned int *ranges;
  unsigned trellis_size;
  uint64_t *trellis_profile;
  void init(unsigned int n, unsigned int k, int **checkmatrix);
 public:
  TrellisDecoder(unsigned int n, unsigned int k, std::ifstream &filename);
  TrellisDecoder(unsigned int n, unsigned int k, int **checkmatrix);
  double decode(const double *x, int *u) override = 0;
  ~TrellisDecoder() override;
};
#endif //BEAST_TRELLISDECODER_H
