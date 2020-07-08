#ifndef BEAST_INCLUDE_KTKLDECODER_H_
#define BEAST_INCLUDE_KTKLDECODER_H_

#include <utility>
#include <bitset>

#include "SoftDecoder.h"
#include "ViterbiDecoder.h"
#include "OrderedStatisticsDecoder.h"

class KTKLDecoder : public SoftDecoder {
 private:
  ViterbiDecoder viterbi;
  OrderedStatisticsDecoder first_candidate;
  unsigned int order, buf_size;
  int **check, *c_best, *l_weights, **best_words, **l_words, *c_last, *c_prev;
  double *gen_metrics, *best_metrics;
  std::vector<int> weight_profile;
  std::vector<int> sorted_beta_ind;
  double l(int **words, int h);
  void index_set_union(std::bitset<128> &first, const std::bitset<128> &second, int range);
  double record_result(const std::bitset<128> &set, int range);
 public:
  KTKLDecoder(unsigned int n, unsigned int k, int **g, int **h, unsigned int w, unsigned int buf_size);
  ~KTKLDecoder();
  double decode(const double *y, int* u) override;
  unsigned int op_iters;
};
#endif //BEAST_INCLUDE_KTKLDECODER_H_
