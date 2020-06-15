#include <ViterbiDecoder.h>
#include <cmath>

#define HUGE_VAL (__builtin_huge_val ())

ViterbiDecoder::ViterbiDecoder(unsigned int n, unsigned int k, int **g) : SoftDecoder(n, k) {
  trellis.construct_from_gen_matrix(n, k, g);
  prev_layer = new double[trellis.get_max_layer_size()];
  cur_layer = new double[trellis.get_max_layer_size()];
  prev_layer_path = new uint64_t[trellis.get_max_layer_size()];
  cur_layer_path = new uint64_t[trellis.get_max_layer_size()];
}
double ViterbiDecoder::decode(double *x, int *u) {
  op_add = 0;
  op_cmp = 0;

  for (unsigned int i = 0; i < n; ++i) {
    alpha[i] = x[i] < 0 ? 0 : 1;
    ++op_cmp;
    beta[i] = fabs(x[i]);
  }

  prev_layer[0] = 0;
  prev_layer_path[0] = 0;
  for (unsigned j = 0; j < n; j++) {
    for (unsigned l = 0; l < trellis[j+1].size(); l++) {
      cur_layer[l] = HUGE_VAL;
      cur_layer_path[l] = 0;
    }
    for (unsigned l = 0; l < trellis[j].size(); l++) {
      for (unsigned z = 0; z < 2; z++) {
        unsigned long long s1 = trellis[j][l].next_node[z];
        if (s1 != ~0) {
          double metric;
          if (z ^ alpha[j]) {
            metric = prev_layer[l] + beta[j];
            ++op_add;
          } else {
            metric = prev_layer[l];
          }
          if (metric < cur_layer[s1]) {
            cur_layer[s1] = metric;
            cur_layer_path[s1] = prev_layer_path[l] ^ (uint64_t(z) << (j));
          }
          ++op_cmp;
        };
      };
    }
    std::swap(prev_layer, cur_layer);
    std::swap(prev_layer_path, cur_layer_path);
  };
  for (unsigned int i = 0; i < n; ++i) {
    u[i] = (prev_layer_path[0] & (uint64_t(1) << i)) ? 1 : 0;
  }
  return prev_layer[0];
}
ViterbiDecoder::~ViterbiDecoder() {

  delete[] prev_layer;
  delete[] cur_layer;
  delete[] prev_layer_path;
  delete[] cur_layer_path;
}
