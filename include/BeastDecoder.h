#ifndef BEAST_BEASTDECODER_H
#define BEAST_BEASTDECODER_H

#include <fstream>
#include <cstdint>
#include <set>
#include <queue>

#include "TrellisDecoder.h"

enum InsertionStatus {
  INSERTED,
  REPLACED,
  DISCARDED,
  ENDED,
};
class BeastDecoder : public TrellisDecoder {
 protected:
  double min_metric = -1;
  double delta = 1;
  uint64_t min_candidate = 0;
  Node *fwd_tree{}, *fwd_tree_buffer{}, *bkw_tree{}, *bkw_tree_buffer{};
  unsigned fwd_tree_size{}, fwd_tree_buffer_size{}, bkw_tree_size{}, bkw_tree_buffer_size{};
  Node **trellis{};

  InsertionStatus insert_node(const Node &node);
  void init(double delta);
 public:
  BeastDecoder(unsigned int n, unsigned int k, std::ifstream &filename, double delta);
  BeastDecoder(unsigned int n, unsigned int k, int **h, double delta);
  double decode(double *x, int *u) override;
  void set_delta(double d);
  double get_delta();
  ~BeastDecoder();
};

#endif //BEAST_BEASTDECODER_H
