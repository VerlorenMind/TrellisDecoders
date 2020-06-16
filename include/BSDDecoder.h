#ifndef BEAST_BSDDECODER_H
#define BEAST_BSDDECODER_H

#include "TrellisDecoder.h"

class BSDDecoder : public TrellisDecoder {
 protected:
  Node *fwd_tree;
  Node *bkw_tree;

  void bubble_insert(Node *array, unsigned start, unsigned end, Node elem);
 public:
  BSDDecoder(unsigned int n, unsigned int k, std::ifstream &filename);
  BSDDecoder(unsigned int n, unsigned int k, int **h);
  ~BSDDecoder();
  double decode(const double *x, int *u) override;
};
#endif //BEAST_BSDDECODER_H
