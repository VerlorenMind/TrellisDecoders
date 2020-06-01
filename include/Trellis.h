#ifndef BEAST_INCLUDE_TRELLIS_H_
#define BEAST_INCLUDE_TRELLIS_H_

#include <cstdint>

struct TrellisNode {
  uint64_t next_node[2];
  uint64_t prev_node[2];
};

class TrellisLayer {
 private:
  TrellisNode* layer;
  unsigned long layer_size;
  unsigned int layer_num;
 public:
  ~TrellisLayer();
  void init(unsigned int size, unsigned int num);
  unsigned int size();
  TrellisNode& operator[](unsigned int i);
};

class Trellis {
 private:
  TrellisLayer *trellis;
  unsigned long long max_layer_size;
 public:
  Trellis();
  ~Trellis();
  // void construct_from_check_matrix(unsigned int n, unsigned int k, int** h);
  void construct_from_gen_matrix(unsigned int n, unsigned int k, int**g);
  TrellisLayer& operator[](unsigned int i);
  unsigned long long get_max_layer_size();
};
#endif //BEAST_INCLUDE_TRELLIS_H_
