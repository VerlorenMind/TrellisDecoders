#include <KTKLDecoder.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iterator>
#include <set>
#include <unordered_set>
#include "CatchWrap.h"

const double EPSILON = 0.000000001;

KTKLDecoder::KTKLDecoder(unsigned int n, unsigned int k, int **g, int **h,
                         unsigned int order, unsigned int buf_size)
    : buf_size(buf_size),
      order(order),
      first_candidate(n, k, g, 1),
      viterbi(n, k, h),
      SoftDecoder(n, k) {
  id = DecoderID::KTKL;
  if (buf_size > 3) buf_size = 3;
  if (buf_size < 1) buf_size = 1;
  gen_metrics = new double[3];
  best_words = new int *[buf_size];
  best_metrics = new double[buf_size];
  l_weights = new int[buf_size];
  l_words = new int *[buf_size];
  c_best = new int[n];
  c_last = new int[n];
  c_prev = new int[n];
  for (unsigned int i = 0; i < buf_size; ++i) {
    best_words[i] = new int[n];
  }
  check = new int *[n - k];
  for (unsigned int i = 0; i < n - k; ++i) {
    check[i] = new int[n];
    memcpy(check[i], h[i], sizeof(int) * n);
  }
  w_profile = weight_profile(n, k, g);
  viterbi.trellis.reduce_to_weight(w_profile[order]);
}
KTKLDecoder::~KTKLDecoder() {
  delete[] gen_metrics;
  for (unsigned int i = 0; i < buf_size; ++i) {
    delete[] best_words[i];
  }
  delete[] best_words;
  delete[] best_metrics;
  delete[] c_last;
  delete[] c_prev;
  for (unsigned int i = 0; i < n - k; ++i) {
    delete[] check[i];
  }
  delete[] check;
  delete[] l_weights;
  delete[] l_words;
  delete[] c_best;
}
#ifdef CATCH_TESTING
bool optimality_holds = false;
#endif

double KTKLDecoder::decode(const double *y, int *u) {
  op_add = 0;
  op_cmp = 0;
  op_iters = 0;
  memset(gen_metrics, 0, 3 * sizeof(double));

#ifdef CATCH_TESTING
  optimality_holds = false;
#endif
  for (unsigned int i = 0; i < n; ++i) {
    alpha[i] = y[i] < 0 ? 0 : 1;
    ++op_cmp;
    beta[i] = fabs(y[i]);
  }

  // Calculating the syndrome of the hard decision
  unsigned int synd_weight = 0;
  for (unsigned int j = 0; j < n - k; ++j) {
    int temp = 0;
    for (unsigned int i = 0; i < n; ++i) {
      temp ^= alpha[i] & check[j][i];
    }
    synd_weight += temp;
  }
  if (synd_weight == 0) {
    for (unsigned int i = 0; i < n; ++i) {
      u[i] = alpha[i];
    }
    return 0;
  }
  // Sorting betas and storing permutation to use it to calculate l
  sorted_beta_ind.resize(n);
  for (unsigned int i = 0; i < n; ++i) {
    sorted_beta_ind[i] = i;
  }
  std::sort(sorted_beta_ind.begin(), sorted_beta_ind.end(),
            [this](const int &a, const int &b) {
              ++op_cmp;
              return beta[a] < beta[b];
            });
  // Generating initial codeword
  int g = 0;
  int last_metric = 0;
  double best_metric;
  best_metric = gen_metrics[last_metric] = first_candidate.decode(y, c_last);
  memcpy(c_best, c_last, n * sizeof(int));
  op_cmp += first_candidate.op_cmp;
  op_add += first_candidate.op_add;
  unsigned int h = 1;
  unsigned int best_words_count = 0;
  unsigned int worst_best_word = 0;
  bool replaced_best = true;
  while (true) {
    ++op_iters;
    // Test optimality of c_best
    double l1;
    for (unsigned int i = 0; i < h - 1; ++i) {
      if (i == h-2 && !replaced_best) {
        l_words[i] = c_best;
      } else {
        l_words[i] = best_words[i];
      }
      l_weights[i] = w_profile[order + 1];
    }
    l_words[h - 1] = c_last;
    l_weights[h - 1] = w_profile[1];
    l1 = l(l_words, h);
    double l2;
    if (g >= buf_size) {
      for (unsigned int i = 0; i < h; ++i) {
        if (i == h-1 && !replaced_best) {
          l_words[i] = c_best;
        } else {
          l_words[i] = best_words[i];
        }
        l_weights[i] = w_profile[order + 1];
      }
      l2 = l(l_words, h);
      l1 = l1 > l2 ? l1 : l2;
    }
    ++op_cmp;
    if (best_metric <= l1) {
#ifdef CATCH_TESTING
      optimality_holds = true;
#endif
      break;
    }
    ++op_cmp;
    // Updating bG
    if (!replaced_best) {
      // If there is place in the buffer, insert c_last in the end
      if(best_words_count < buf_size) {
        memcpy(best_words[best_words_count], c_last, n * sizeof(int));
        best_metrics[best_words_count] = gen_metrics[last_metric];
        if(best_metrics[best_words_count] > best_metrics[worst_best_word]) worst_best_word = best_words_count;
        ++best_words_count;
      }
      // Else, checking if it is better than the worst best word (what)
      else if(best_metrics[worst_best_word] >= gen_metrics[last_metric]) {
        for(unsigned int j = worst_best_word; j < best_words_count-1; ++j) {
          memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
          best_metrics[j] = best_metrics[j + 1];
        }
        memcpy(best_words[best_words_count], c_last, n * sizeof(int));
        best_metrics[best_words_count] = gen_metrics[last_metric];
        worst_best_word = 0;
        for(unsigned int j=1; j<best_words_count; ++j) {
          if(best_metrics[worst_best_word] < best_metrics[j]) {
            worst_best_word = j;
          }
        }
      }
    }
    // If not optimal, search trellis
    // Searching subtrellis
    if (last_metric < 2) {
      ++last_metric;
    } else {
      memmove(gen_metrics, gen_metrics + 1, 2 * sizeof(double));
    }
    std::swap(c_last, c_prev);
    ++g;
    h = g + 1 < buf_size ? g + 1 : buf_size;
    gen_metrics[last_metric] = viterbi.decode_around(y, c_last, c_prev);
    op_cmp += viterbi.op_cmp;
    op_add += viterbi.op_add;
    // Updating c_best
    replaced_best = false;
    if (gen_metrics[last_metric] < best_metric) {
      if(best_words_count >= buf_size) {
        for (int j = 0; j < best_words_count-1; ++j) {
          memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
          best_metrics[j] = best_metrics[j + 1];
        }
        worst_best_word = 1;
        for(unsigned int j=2; j<best_words_count; ++j) {
          if(best_metrics[worst_best_word] < best_metrics[j]) {
            worst_best_word = j;
          }
        }
      }
      memcpy(best_words[best_words_count], c_best, n * sizeof(int));
      best_metrics[best_words_count] = best_metric;
      best_words_count = best_words_count == buf_size ? best_words_count : best_words_count + 1;
      memcpy(c_best, c_last, n * sizeof(int));
      best_metric = gen_metrics[last_metric];
      replaced_best = true;
    }
    // Test stopping condition
    if (g >= 2) {
      if (gen_metrics[2] - gen_metrics[0] < EPSILON) {
        break;
      }
    }
  }
  for (unsigned int i = 0; i < n; ++i) {
    u[i] = c_best[i];
  }
  return best_metric;
}
double KTKLDecoder::l(int **words, int h) {
  // I don't know what's happening either
  double res = HUGE_VAL;
  if (h == 1) {
    std::bitset<128> d0, d1;
    d0.reset();
    d1.reset();
    for (unsigned int i = 0; i < n; ++i) {
      if (words[0][i] == alpha[i]) {
        d0.set(i);
      } else {
        d1.set(i);
      }
    }
    int delta = l_weights[0] - (int) d1.count();
    res = record_result(d0, delta);
  } else if (h == 2) {
    std::bitset<128> d0[2], d1[2];
    for (unsigned int z = 0; z < 2; ++z) {
      d0[z].reset();
      d1[z].reset();
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] == alpha[i]) {
          d0[z].set(i);
        } else {
          d1[z].set(i);
        }
      }
    }
    int delta1 = l_weights[0] - (int) d1[0].count();
    int delta2 = l_weights[1] - (int) d1[1].count();
    if (delta1 < delta2) {
      std::swap(delta1, delta2);
      std::swap(d0[0], d0[1]);
      std::swap(d1[0], d1[1]);
    }
    std::bitset<128> d00 = d0[0] & d0[1];
    std::bitset<128> d01 = d0[0] & d1[1];
    int range = (delta1 - delta2) / 2;
    index_set_union(d00, d01, range);
    res = record_result(d00, delta1);
  } else if (h == 3) {
    double l1, l2;
    std::bitset<128> d0[3], d1[3];
    std::bitset<128> d[8];
    int delta[3];
    for (unsigned int z = 0; z < 3; ++z) {
      d0[z].reset();
      d1[z].reset();
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] == alpha[i]) {
          d0[z].set(i);
        } else {
          d1[z].set(i);
        }
      }
      delta[z] = l_weights[z] - (int) d1[z].count();
    }
    std::vector<int> delta_order;
    delta_order.resize(3);
    for (int i = 0; i < 3; ++i) delta_order[i] = i;
    std::sort(delta_order.begin(), delta_order.end(), [&delta](const int &a, const int &b) {
      return delta[a] > delta[b];
    });
    for (unsigned int i = 0; i < 3; ++i) {
      while (delta_order[i] != i) {
        std::swap(d0[i], d0[delta_order[i]]);
        std::swap(d1[i], d1[delta_order[i]]);
        std::swap(delta[i], delta[delta_order[i]]);
        std::swap(delta_order[i], delta_order[delta_order[i]]);
      }
    }
    for (unsigned int i = 0; i < 8; ++i) {
      d[i].reset();
      if (i & 1) {
        d[i] = d1[0];
      } else {
        d[i] = d0[0];
      }
      if (i & 2) {
        d[i] &= d1[1];
      } else {
        d[i] &= d0[1];
      }
      if (i & 4) {
        d[i] &= d1[2];
      } else {
        d[i] &= d0[2];
      }
    }
    int delta12 = (delta[0] - delta[1]) / 2;
    delta12 = std::min(delta[0], delta12);
    int delta13 = (delta[0] - delta[2]) / 2;
    delta13 = std::min(delta[0], delta13);
    int delta1 = 0;
    delta1 = std::max(delta12, (int) d[0b010].count());
    delta1 = std::max(delta13 - (int) d[0b001].count(), delta1);
    int delta2 = (int) d[0b011].count();
    delta2 = std::min(delta12, delta2);
    delta2 = std::min(delta13, delta2);
    delta2 = std::min((int) d[0b000].count() + delta12 + delta13 - delta[0], delta2);
    l1 = HUGE_VAL;
    for (int i = delta1; i <= delta2; ++i) {
      std::bitset<128> I, I2;
      I.reset();
      I2.reset();
      if (delta[0] - i > 0) {
        I2 = d[0b000];
        index_set_union(I2, d[0b001], delta13 - i);
        index_set_union(I2, d[0b010], delta12 - i);
        index_set_union(I, I2, delta[0] - i);
      }
      index_set_union(I, d[0b011], i);
      double temp = record_result(I, (int) I.count());
      if (l1 > temp) {
        l1 = temp;
      }
      ++op_cmp;
    }
    int eps[3] = {0, 0, 0};
    if (!(delta[0] % 2 == delta[1] % 2 && delta[1] % 2 == delta[2] % 2)) {
      if (delta[0] % 2 != delta[1] % 2 && delta[0] % 2 != delta[2]) {
        eps[0] = 1;
      } else if (delta[1] % 2 != delta[0] % 2) {
        eps[1] = 1;
      } else {
        eps[2] = 1;
      }
    }
    delta12 = (delta[0] + eps[0] - delta[1] - eps[1]) / 2;
    delta13 = (delta[0] + eps[0] - delta[2] - eps[2]) / 2;
    int mean_delta12 =
        ((delta[1] + delta[2]) / 2) +
            ((delta[1] + delta[2]) % 2);  // ceil((delta[1] + delta[2]) / 2)
    delta1 = mean_delta12 - (int) d[0b000].count();
    delta1 = std::max(delta1, 0);
    delta2 = (int) d[0b100].count();
    delta2 = std::min(delta2, ((int) d[0b010].count() - delta12));
    delta2 = std::min(delta2, ((int) d[0b001].count() - delta13));
    delta2 = std::min(delta2, mean_delta12);
    l2 = HUGE_VAL;
    for (int i = delta1; i <= delta2; ++i) {
      std::bitset<128> I;
      index_set_union(I, d[0b000], mean_delta12 - i);
      index_set_union(I, d[0b010], delta12 + i);
      index_set_union(I, d[0b001], delta13 + i);
      index_set_union(I, d[0b100], i);
      double temp = record_result(I, (int) I.count());
      if (l2 > temp) {
        l2 = temp;
      }
      ++op_cmp;
    }
    res = std::min(l1, l2);
    ++op_cmp;
  } else {
    res = HUGE_VAL;
  }
  return res;
}

void KTKLDecoder::index_set_union(std::bitset<128> &first,
                                  const std::bitset<128> &second,
                                  int range) {
  auto iter = sorted_beta_ind.begin();
  for (int i = 0; i < second.count() && i < range; ++i) {
    for (; iter != sorted_beta_ind.end() && !second[*iter]; ++iter) {
    }
    assert(iter != sorted_beta_ind.end());
    first.set(*iter);
    ++iter;
  }
}

double KTKLDecoder::record_result(const std::bitset<128> &set, int range) {
  auto iter = sorted_beta_ind.begin();
  double temp = 0;
  for (int i = 0; i < set.count() && i < range; ++i) {
    for (; iter != sorted_beta_ind.end() && !set[*iter]; ++iter) {
    }
    assert(iter != sorted_beta_ind.end());
    temp += beta[*iter];
    ++op_add;
    ++iter;
  }
  return temp;
}
