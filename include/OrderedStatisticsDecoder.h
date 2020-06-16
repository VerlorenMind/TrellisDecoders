#ifndef BEAST_INCLUDE_ORDEREDSTATISTICSDECODER_H_
#define BEAST_INCLUDE_ORDEREDSTATISTICSDECODER_H_

#include "SoftDecoder.h"

class OrderedStatisticsDecoder : public SoftDecoder {
 private:
  int **g;
  int **g_buf;
  int *c_temp;
  int *c_best;
  int *u_temp;
  int w;
 public:
  OrderedStatisticsDecoder(unsigned int n, unsigned int k, int **gen_matrix, int w);
  ~OrderedStatisticsDecoder();
  double decode(const double *x, int *u) override;
  void set_max_weight(int w);
};
#endif //BEAST_INCLUDE_ORDEREDSTATISTICSDECODER_H_
