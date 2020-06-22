#ifndef BEAST_INCLUDE_TRELLIS_H_
#define BEAST_INCLUDE_TRELLIS_H_

#include <cstdint>
#include <ostream>
#include <vector>

struct TrellisNode {
  uint64_t next_node[2];
  uint64_t prev_node[2];
};

class TrellisLayer {
 private:
  TrellisNode* layer;
  unsigned long layer_size;
  unsigned long true_size;
 public:
  ~TrellisLayer();
  void init(unsigned int size, unsigned int num);
  unsigned int size();
  TrellisNode& operator[](unsigned int i);
  void delete_node(unsigned long long num);
  void add_node();
  void shrink_to_fit();
};

class Trellis {
 private:
  TrellisLayer *trellis;
  unsigned long long max_layer_size;
  unsigned long trellis_size;

  bool dfs_trellis_search(std::pair<unsigned int, unsigned long long> node,
                          std::vector<bool> *keep_node, unsigned int metric, unsigned int w,
                          std::vector<int> &path, std::vector<std::vector<int>> &paths_to_keep);
  unsigned long long second_0_edge_to_end = ~0;
 public:
  Trellis();
  ~Trellis();
  void construct_from_check_matrix(unsigned int n, unsigned int k, int** h);
  void construct_from_gen_matrix(unsigned int n, unsigned int k, int**g);
  TrellisLayer& operator[](unsigned int i);
  unsigned long long get_max_layer_size();
  void print_trellis(std::ostream &out);
  void reduce_to_weight(unsigned int w);
  void reduce_to_weight_dfs(unsigned int w);
  void delete_node(unsigned int layer, unsigned long long num);
};
#endif //BEAST_INCLUDE_TRELLIS_H_
