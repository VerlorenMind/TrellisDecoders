#include <KTKLDecoder.h>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iterator>
#include <set>
#include <unordered_set>
#include "CatchWrap.h"

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
  weight_profile.resize(0);
  int word_weight = 0;
  int *codeword = new int[n];
  for (uint64_t infoword = 0; infoword < (uint64_t(1) << k); ++infoword) {
    if (infoword == 0) {
      memset(codeword, 0, n * sizeof(int));
    } else {
      uint64_t bit_diff = infoword ^ (infoword - 1);
      for (unsigned int l = 0; l < k; ++l) {
        if (bit_diff & (uint64_t(1) << l)) {
          for (unsigned int j = 0; j < n; ++j) {
            codeword[j] ^= g[l][j];
          }
        }
      }
    }
    word_weight = 0;
    for (unsigned int i = 0; i < n; ++i) {
      word_weight += codeword[i];
    }
    auto iter = weight_profile.begin();
    while (iter != weight_profile.end() && *iter < word_weight) {
      ++iter;
    }
    if (iter == weight_profile.end()) {
      weight_profile.push_back(word_weight);
    } else if (*iter != word_weight) {
      weight_profile.insert(iter, word_weight);
    }
  }
  delete[] codeword;
  viterbi.trellis.reduce_to_weight(weight_profile[order]);
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
              return beta[a] <= beta[b];
            });
  // Generating initial codeword
  int g = -1;
  int last_metric = 0;
  double best_metric;
  best_metric = gen_metrics[last_metric] = first_candidate.decode(y, c_last);
  memcpy(c_best, c_last, n * sizeof(int));
  op_cmp += first_candidate.op_cmp;
  op_add += first_candidate.op_add;
  unsigned int h = 1;
  unsigned int best_words_count = 0;
  bool replaced_best = true;
  while (true) {
    // Test optimality of c_best
    double l1;
    for (unsigned int i = 0; i < h - 1; ++i) {
      l_words[i] = best_words[i];
      l_weights[i] = weight_profile[order + 1];
    }
    l_words[h - 1] = c_last;
    l_weights[h - 1] = weight_profile[1];
    l1 = l(l_words, h);
    double l2;
    if (g >= buf_size) {
      for (unsigned int i = 0; i < h; ++i) {
        l_words[i] = best_words[i];
        l_weights[i] = weight_profile[order + 1];
      }
      l2 = l(l_words, h);
    } else {
      l2 = 0;
    }
    l1 = l1 > l2 ? l1 : l2;
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
      // If buffer is empty, just put the word there
      if (best_words_count == 0) {
        memcpy(best_words[0], c_last, n * sizeof(int));
        best_metrics[0] = gen_metrics[last_metric];
        ++best_words_count;
      } else {
        // Finding the place to insert the word
        unsigned int i = 0;
        for (i = 0;
             i < best_words_count && best_metrics[i] > gen_metrics[last_metric];
             ++i) {
          ++op_cmp;
        }
        ++op_cmp;
        // If there is still space in the buffer, make room for the new word
        if (best_words_count < buf_size) {
          for (int j = (int)best_words_count; j > (int)i; --j) {
            memcpy(best_words[j], best_words[j - 1], n * sizeof(int));
            best_metrics[j] = best_metrics[j - 1];
          }
          ++best_words_count;
          memcpy(best_words[i], c_last, n * sizeof(int));
          best_metrics[i] = gen_metrics[last_metric];
        }
        // If not, if there's inferior words in the buffer with higher metric,
        // delete the one with the highest and insert the new one Otherwise,
        // don't insert anything
        else if (i > 0) {
          for (unsigned int j = 0; j < i - 1; ++j) {
            memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
            best_metrics[j] = best_metrics[j + 1];
          }
          memcpy(best_words[i - 1], c_last, n * sizeof(int));
          best_metrics[i - 1] = gen_metrics[last_metric];
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
    h = g+1 < buf_size ? g+1 : buf_size;
    gen_metrics[last_metric] = viterbi.decode_around(y, c_last, c_prev);
    op_cmp += viterbi.op_cmp;
    op_add += viterbi.op_add;
    // Updating c_best
    replaced_best = false;
    if (gen_metrics[last_metric] < best_metric) {
      if (best_words_count == buf_size) {
        for (unsigned int j = 0; j < best_words_count - 1; ++j) {
          memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
          best_metrics[j] = best_metrics[j + 1];
        }
      } else {
        ++best_words_count;
      }
      memcpy(best_words[best_words_count - 1], c_best, n * sizeof(int));
      best_metrics[best_words_count - 1] = best_metric;
      memcpy(c_best, c_last, n * sizeof(int));
      best_metric = gen_metrics[last_metric];
      replaced_best = true;
    }
    ++op_cmp;
    // Test stopping condition
    if (g >= 2) {
      if (gen_metrics[2] == gen_metrics[0]) {
        break;
      }
      ++op_cmp;
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
    for (unsigned int i = 0; i < n; ++i) {
      if (words[0][i] == alpha[i]) {
        d0.set(i);
      } else {
        d1.set(i);
      }
    }
    int delta = l_weights[0] - (int)d1.count();
    res = 0;
    res = record_result(d0, delta);
  } else if (h == 2) {
    std::bitset<128> d0[2], d1[2];
    d0[0].reset();
    d0[1].reset();
    d1[0].reset();
    d1[1].reset();
    for (unsigned int z = 0; z < 2; ++z) {
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] != alpha[i]) {
          d0[z].set(i);
        } else {
          d1[z].set(1);
        }
      }
    }
    int delta1 = l_weights[0] - (int)d1[0].count();
    int delta2 = l_weights[1] - (int)d1[1].count();
    if (delta1 < delta2) {
      std::swap(delta1, delta2);
      std::swap(d0[0], d0[1]);
      std::swap(d1[0], d1[1]);
    }
    std::bitset<128> d00, d01;
    d00 = d0[0] & d0[1];
    d01 = d0[0] & d1[1];
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
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] != alpha[i]) {
          d0[z].set(i);
        } else {
          d1[z].set(i);
        }
      }
      delta[z] = l_weights[z] - (int)d1[z].count();
    }
    std::vector<int> order;
    order.resize(3);
    for (int i = 0; i < 3; ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&delta](int a, int b) {
      return delta[a] >= delta[b]; });
    for (unsigned int i = 0; i < 3; ++i) {
      while (order[i] != i) {
        std::swap(d0[i], d0[order[i]]);
        std::swap(d1[i], d1[order[i]]);
        std::swap(delta[i], delta[order[i]]);
        std::swap(order[i], order[order[i]]);
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
    delta12 = delta[0] < delta12 ? delta[0] : delta12;
    int delta13 = (delta[0] - delta[2]) / 2;
    delta13 = delta[0] < delta13 ? delta[0] : delta13;
    int delta1 = 0;
    delta1 = delta12 - (int)d[0b010].count() > delta1 ? delta12 - (int)d[0b010].count()
                                                     : delta1;
    delta1 = delta13 - (int)d[0b001].count() > delta1 ? delta13 - (int)d[0b001].count()
                                                     : delta1;
    int delta2 = (int)d[0b011].count();
    delta2 = delta12 < delta2 ? delta12 : delta2;
    delta2 = delta13 < delta2 ? delta13 : delta2;
    delta2 = (int)d[0b000].count() + delta12 + delta13 - delta[0] < delta2
                 ? (int)d[0b000].count() + delta12 + delta13 - delta[0]
                 : delta2;
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
      double temp = record_result(I, (int)I.count());
      if (l1 < temp) {
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
    delta1 = mean_delta12;
    delta1 -= (int)d[0b000].count();
    delta1 = delta1 > 0 ? delta1 : 0;
    delta2 = (int)d[0b100].count();
    delta2 = delta2 < ((int)d[0b010].count() - delta12)
                 ? delta2
                 : (int)d[0b010].size() - delta12;
    delta2 = delta2 < ((int)d[0b001].count() - delta13)
                 ? delta2
                 : (int)d[0b001].count() - delta13;
    delta2 = delta2 < mean_delta12 ? delta2 : mean_delta12;
    l2 = HUGE_VAL;
    for (int i = delta1; i <= delta2; ++i) {
      std::bitset<128> I;
      index_set_union(I, d[0b000], mean_delta12 - i);
      index_set_union(I, d[0b010], delta12 + i);
      index_set_union(I, d[0b001], delta13 + i);
      index_set_union(I, d[0b100], i);
      double temp = record_result(I, (int)I.count());
      if (l1 < temp) {
        l1 = temp;
      }
      ++op_cmp;
    }
    res = l1 < l2 ? l1 : l2;
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
