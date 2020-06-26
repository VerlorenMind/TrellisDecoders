
#include <KTKLDecoder.h>
#include <cmath>
#include <set>
#include <algorithm>
#include <unordered_set>

KTKLDecoder::KTKLDecoder(unsigned int n, unsigned int k, int **g, int **h, unsigned int order, unsigned int buf_size)
    : buf_size(buf_size), order(order), first_candidate(n, k, g, 1), viterbi(n, k, h), SoftDecoder(n, k) {
  id = DecoderID::KTKL;
  if (buf_size > 3) buf_size = 3;
  if (buf_size < 1) buf_size = 1;
  gen_metrics = new double[buf_size + 1];
  gen_words = new int *[buf_size + 1];
  best_words = new int *[buf_size];
  best_metrics = new double[buf_size];
  l_weights = new int[buf_size];
  l_words = new int*[buf_size];
  c_best = new int[n];
  for (unsigned int i = 0; i < buf_size; ++i) {
    gen_words[i] = new int[n];
    best_words[i] = new int[n];
  }
  gen_words[buf_size] = new int[n];
  check = new int *[n - k];
  for (unsigned int i = 0; i < n - k; ++i) {
    check[i] = new int[n];
    memcpy(check[i], h[i], sizeof(int) * n);
  }
  weight_profile.resize(0);
  int word_weight = 0;
  for (uint64_t infoword = 0; infoword < (uint64_t(1) << k); ++infoword) {
    word_weight = 0;
    for (unsigned int j = 0; j < n; ++j) {
      int temp = 0;
      for (unsigned int l = 0; l < k; ++l) {
        temp ^= g[l][j] & ((infoword & (uint64_t(1) << l)) ? 1 : 0);
      }
      word_weight += temp ? 1 : 0;
    }
    auto iter = weight_profile.begin();
    while (iter != weight_profile.end() && *iter < word_weight) {
      ++iter;
    }
    if(iter == weight_profile.end()) {
      weight_profile.push_back(word_weight);
    }
    else if (*iter != word_weight) {
      weight_profile.insert(iter, word_weight);
    }
  }
  viterbi.trellis.reduce_to_weight(weight_profile[order]);
}
KTKLDecoder::~KTKLDecoder() {
  delete[] gen_metrics;
  for (unsigned int i = 0; i < buf_size; ++i) {
    delete[] gen_words[i];
    delete[] best_words[i];
  }
  delete[] best_words;
  delete[] best_metrics;
  delete[] gen_words[buf_size];
  delete[] gen_words;
  for (unsigned int i = 0; i < n - k; ++i) {
    delete[] check[i];
  }
  delete[] check;
  delete[] l_weights;
  delete[] l_words;
  delete[] c_best;
}
double KTKLDecoder::decode(const double *y, int *u) {
  op_add = 0;
  op_cmp = 0;

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
              return beta[a] > beta[b];
            });
  // Generating initial codeword
  int g = 0;
  double best_metric;
  best_metric = gen_metrics[g] = first_candidate.decode(y, gen_words[g]);
  memcpy(c_best, gen_words[g], n * sizeof(int));
  op_cmp += first_candidate.op_cmp;
  op_add += first_candidate.op_add;
  int best_words_count = 0;
  bool replaced_best = true;
  while (true) {
    // Test optimality of c_best
    double l1;
    if(best_words_count != 0) {
      for (unsigned int i = 0; i < best_words_count - 1; ++i) {
        l_words[i] = best_words[i];
        l_weights[i] = weight_profile[order + 1];
      }
      l_words[best_words_count - 1] = gen_words[g];
      l_weights[best_words_count - 1] = weight_profile[1];
      l1 = l(l_words, best_words_count);
    }
    else {
      l_words[0] = gen_words[g];
      l_weights[0] = weight_profile[1];
      l1 = l(l_words, best_words_count+1);
    }
    for(unsigned int i = 0; i < best_words_count; ++i) {
      l_words[i] = best_words[i];
      l_weights[i] = weight_profile[order+1];
    }
    double l2;
    if(g == buf_size) {
      l2 = l(l_words, best_words_count);
    }
    else {
      l2 = 0;
    }
    l1 = l1 > l2 ? l1 : l2;
    ++op_cmp;
    if(best_metric < l1) {
      break;
    }
    ++op_cmp;
    // Updating bG
    if (!replaced_best) {
      // If buffer is empty, just put the word there
      if (best_words_count == 0) {
        memcpy(best_words[0], gen_words[g], n * sizeof(int));
        best_metrics[0] = gen_metrics[g];
        ++best_words_count;
      } else {
        // Finding the place to insert the word
        unsigned int i = 0;
        for (i = 0; i < best_words_count && best_metrics[i] > gen_metrics[g]; ++i) {
          ++op_cmp;
        }
        ++op_cmp;
        // If there is still space in the buffer, make room for the new word
        if (best_words_count < buf_size) {
          for (int j = best_words_count; j > i; --j) {
            memcpy(best_words[j], best_words[j - 1], n * sizeof(int));
            best_metrics[j] = best_metrics[j - 1];
          }
          ++best_words_count;
          memcpy(best_words[i], gen_words[g], n * sizeof(int));
          best_metrics[i] = gen_metrics[g];
        }
        // If not, if there's inferior words in the buffer with higher metric, delete the one with the highest and insert the new one
        // Otherwise, don't insert anything
        else if (i > 0) {
          for (int j = 0; j < i - 1; ++j) {
            memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
            best_metrics[j] = best_metrics[j + 1];
          }
          memcpy(best_words[i - 1], gen_words[g], n * sizeof(int));
          best_metrics[i - 1] = gen_metrics[g];
        }
      }
    }
    // If not optimal, search trellis
    // If buffer is full, shifting saved words to make space for a new one
    if (g == buf_size) {
      for (int i = buf_size; i > 0; --i) {
        memcpy(gen_words[i - 1], gen_words[i], n * sizeof(int));
      }
    } else {
      ++g;
    }
    // Searching subtrellis
    gen_metrics[g] = viterbi.decode_around(y, gen_words[g], gen_words[g-1]);
    op_cmp += viterbi.op_cmp;
    op_add += viterbi.op_add;
    // Test stopping condition
    if (g >= 2) {
      if (gen_metrics[g] == gen_metrics[g - 2]) {
        break;
      }
      ++op_cmp;
    }
    // Updating c_best
    replaced_best = false;
    if (gen_metrics[g] < best_metric) {
      if (best_words_count == buf_size) {
        for (int j = 0; j < best_words_count - 1; ++j) {
          memcpy(best_words[j], best_words[j + 1], n * sizeof(int));
          best_metrics[j] = best_metrics[j + 1];
        }
      } else {
        ++best_words_count;
      }
      memcpy(best_words[best_words_count - 1], c_best, n * sizeof(int));
      memcpy(c_best, gen_words[g], n * sizeof(int));
      replaced_best = true;
    }
    ++op_cmp;
  }
  for (unsigned int i = 0; i < n; ++i) {
    u[i] = c_best[i];
  }
  return best_metric;
}
double KTKLDecoder::l(int **words, int h) {
  auto ind_comp = [&temp = sorted_beta_ind](const int &a, const int &b) {
    int pos_a = ~0, pos_b = ~0;
    for (unsigned int i = 0; i < temp.size(); ++i) {
      if (temp[i] == a) {
        pos_a = i;
      }
      if (temp[i] == b) {
        pos_b = i;
      }
    }
    return pos_a < pos_b;
  };
  using Set = std::set<int, decltype(ind_comp)>;
  // I don't know what's happening either
  double res = HUGE_VAL;
  if (h == 1) {
    Set d0(ind_comp), d1(ind_comp);
    for (unsigned int i = 0; i < n; ++i) {
      if (words[0][i] != alpha[i]) {
        d1.insert(i);
      } else {
        d0.insert(i);
      }
    }
    int delta = l_weights[0] - d1.size();
    res = 0;
    int range = delta >= d0.size() ? d0.size() : delta;
    for (unsigned int i = 0; i < range; ++i) {
      res += beta[i];
      ++op_add;
    }
  } else if (h == 2) {
    Set d0[2] = {Set(ind_comp), Set(ind_comp)};
    Set d1[2] = {Set(ind_comp), Set(ind_comp)};
    for (unsigned int z = 0; z < 2; ++z) {
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] != alpha[i]) {
          d0[z].insert(i);
        } else {
          d1[z].insert(i);
        }
      }
    }
    int delta1 = l_weights[0] - d1[0].size();
    int delta2 = l_weights[1] - d1[1].size();
    Set d00(ind_comp), d01(ind_comp);
    std::set_intersection(d0[0].begin(), d0[0].end(), d0[1].begin(), d0[1].end(), std::inserter(d00, d00.begin()));
    std::set_intersection(d0[0].begin(), d0[0].end(), d1[1].begin(), d1[1].end(), std::inserter(d01, d01.begin()));
    int range = (delta1 - delta2) / 2 < d01.size() ? (delta1 - delta2) / 2 : d01.size();
    auto iter = d01.begin();
    for (unsigned int i = 0; i < range; ++i) {
      d00.insert(*iter);
      ++iter;
    }
    range = delta1 < d00.size() ? delta1 : d00.size();
    iter = d00.begin();
    for (unsigned int i = 0; i < range; ++i) {
      res += beta[*iter];
      ++op_add;
      ++iter;
    }
  } else if (h == 3) {
    double l1, l2;
    Set d0[3] = {Set(ind_comp), Set(ind_comp), Set(ind_comp)};
    Set d1[3] = {Set(ind_comp), Set(ind_comp), Set(ind_comp)};
    Set d[8] = {Set(ind_comp), Set(ind_comp), Set(ind_comp), Set(ind_comp), Set(ind_comp), Set(ind_comp), Set(ind_comp),
                Set(ind_comp)};
    int delta[3];
    for (unsigned int z = 0; z < 3; ++z) {
      for (unsigned int i = 0; i < n; ++i) {
        if (words[z][i] != alpha[i]) {
          d0[z].insert(i);
        } else {
          d1[z].insert(i);
        }
      }
      delta[z] = l_weights[z] - d1[z].size();
    }
    for (unsigned int i = 0; i < 8; ++i) {
      if (i % 2) {
        std::copy(d1[0].begin(), d1[0].end(), std::inserter(d[i], d[i].begin()));
      } else {
        std::copy(d0[0].begin(), d0[0].end(), std::inserter(d[i], d[i].begin()));
      }
      for (unsigned int j = 1; j < 3; ++j) {
        auto &temp_set = (i & (1 << j)) ? d1[j] : d0[j];
        std::set_intersection(d[i].begin(),
                              d[i].end(),
                              temp_set.begin(),
                              temp_set.end(),
                              std::inserter(d[i], d[i].begin()));
      }
    }
    int delta12 = (delta[0] - delta[1]) / 2;
    delta12 = delta[0] < delta12 ? delta[0] : delta12;
    int delta13 = (delta[0] - delta[2]) / 2;
    delta13 = delta[0] < delta13 ? delta[0] : delta13;
    int delta1 = 0;
    delta1 = delta12 - d[0b010].size() > delta1 ? delta12 - d[0b010].size() : delta1;
    delta1 = delta13 - d[0b001].size() > delta1 ? delta13 - d[0b001].size() : delta1;
    int delta2 = d[0b011].size();
    delta2 = delta12 > delta1 ? delta12 : delta2;
    delta2 = delta13 > delta1 ? delta13 : delta2;
    delta2 = d[0b000].size() + delta12 + delta13 - delta[0] > delta1 ? d[0b000].size() + delta12 + delta13 - delta[0]
                                                                     : delta2;
    if (delta1 <= delta2) {
      l1 = HUGE_VAL;
      for (int i = delta1; i <= delta2; ++i) {
        Set I(ind_comp);
        std::copy(d[0b000].begin(), d[0b000].end(), std::inserter(I, I.begin()));
        int range = (delta13 - i) >= d[0b001].size() ? d[0b001].size() : (delta13 - i);
        auto iter = d[0b001].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        range = (delta12 - i) >= d[0b010].size() ? d[0b010].size() : (delta12 - i);
        iter = d[0b010].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        iter = I.begin();
        if (delta[0] - i < I.size()) {
          for (unsigned int j = 0; j < (delta[0] - j); ++j) {
            ++iter;
          }
        }
        I.erase(iter, I.end());
        range = i < d[0b011].size() ? i : d[0b011].size();
        iter = d[0b011].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        double temp = 0;
        for (int ind : I) {
          temp += beta[ind];
          ++op_add;
        }
        if (l1 < temp) {
          l1 = temp;
        }
        ++op_cmp;
      }
    } else {
      l1 = HUGE_VAL;
    }
    int eps[3];
    eps[0] = eps[1] = eps[2] = 0;
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
    int mean_delta12 = ((delta[1] + delta[2]) / 2) + ((delta[1] + delta[2]) % 2); // ceil((delta[1] + delta[2]) / 2)
    delta1 = mean_delta12;
    delta1 -= d[0b000].size();
    delta1 = delta1 > 0 ? delta1 : 0;
    delta2 = d[0b100].size();
    delta2 = delta2 < (d[0b010].size() - delta12) ? delta2 : d[0b010].size() - delta12;
    delta2 = delta2 < (d[0b001].size() - delta13) ? delta2 : d[0b001].size() - delta13;
    delta2 = delta2 < mean_delta12 ? delta2 : mean_delta12;
    if (delta1 < delta2) {
      l2 = HUGE_VAL;
      for (int i = delta1; i <= delta2; ++i) {
        Set I(ind_comp);
        int range = mean_delta12 - i < d[0b000].size() ? mean_delta12 : d[0b000].size();
        auto iter = d[0b000].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        range = (delta12 + i) < d[0b010].size() ? delta12 + i : d[0b010].size();
        iter = d[0b010].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        range = (delta13 + i) < d[0b001].size() ? delta13 + i : d[0b001].size();
        iter = d[0b001].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        range = i < d[0b100].size() ? i : d[0b100].size();
        iter = d[0b100].begin();
        for (unsigned int j = 0; j < range; ++j) {
          I.insert(*iter);
          ++iter;
        }
        double temp = 0;
        for (int ind : I) {
          temp += beta[ind];
          ++op_add;
        }
        if (l1 < temp) {
          l1 = temp;
        }
        ++op_cmp;
      }
    } else {
      l2 = HUGE_VAL;
    }
    ++op_cmp;
    res = l1 < l2 ? l1 : l2;
  } else {
    res = HUGE_VAL;
  }
  return res;
}

