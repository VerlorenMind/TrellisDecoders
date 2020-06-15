#ifndef BEAST_INCLUDE_KTKLDECODER_H_
#define BEAST_INCLUDE_KTKLDECODER_H_

#include <utility>

#include "SoftDecoder.h"
#include "ViterbiDecoder.h"
#include "OrderedStatisticsDecoder.h"

class KTKLDecoder : public SoftDecoder {
 private:
  ViterbiDecoder viterbi;
  OrderedStatisticsDecoder first_candidate;
  unsigned int order, buf_size;
  int **check, *c_best, *l_weights, **gen_words, **best_words, **l_words;
  double *gen_metrics, *best_metrics, *temp_y;
  std::vector<int> weight_profile;
  std::vector<int> sorted_beta_ind;
  double l(int **words, int h);
 public:
  KTKLDecoder(unsigned int n, unsigned int k, int **g, int **h, unsigned int w, unsigned int buf_size);
  ~KTKLDecoder();
  double decode(double *y, int* u) override;
};
#endif //BEAST_INCLUDE_KTKLDECODER_H_
