#include <OrderedStatisticsDecoder.h>
#include <cstring>
#include <vector>
#include <Utility.h>
#include <cmath>

OrderedStatisticsDecoder::OrderedStatisticsDecoder(unsigned int n, unsigned int k, int **gen_matrix, int w)
    : SoftDecoder(n, k) {
  this->w = w;
  id = DecoderID::ORDERED_STATISTICS;
  g = new int *[k];
  g_buf = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    g[i] = new int[n];
    g_buf[i] = new int[n];
    memcpy(g[i], gen_matrix[i], n * sizeof(int));
  }
  c_temp = new int[n];
  c_best = new int[n];
  c_hard = new int[n];
}

OrderedStatisticsDecoder::~OrderedStatisticsDecoder() {
  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
    delete[] g_buf[i];
  }
  delete[] g;
  delete[] g_buf;
  delete[] c_temp;
  delete[] c_best;
  delete[] c_hard;
}

double OrderedStatisticsDecoder::decode(const double *x, int *u) {
  op_cmp = 0;
  op_add = 0;
  for (unsigned int i = 0; i < n; ++i) {
    beta[i] = fabs(x[i]);
  }
  std::vector<int> permutation;
  op_cmp += quicksort(n, beta, permutation);
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < k; ++j) {
      g_buf[j][i] = g[j][permutation[i]];
    }
  }
  for (unsigned int i = 0; i < n; ++i) {
    alpha[i] = x[permutation[i]] < 0 ? 0 : 1;
    ++op_cmp;
  }
  // Making the matrix diagonal in k first linearly independent columns and storing the indices of these columns
  echelonForm(n, k, g_buf);
  std::vector<int> ind(k);
  int r = 0;
  int c = 0;
  while (r < k && c < n) {
    if (g_buf[r][c]) {
      for (int i = r - 1; i >= 0; --i) {
        if (g_buf[i][c]) {
          for (unsigned int j = 0; j < n; ++j) {
            g_buf[i][j] ^= g_buf[r][j];
          }
        }
      }
      ind[r] = c;
      ++r;
      ++c;
    } else {
      ++c;
    }
  }
  double min_metric = -1;
  memset(c_hard, 0, n*sizeof(int));
  for(unsigned int i=0; i<k; ++i) {
    for(unsigned int j=0; j<n; ++j) {
      c_hard[j] ^= alpha[ind[i]] & g_buf[i][j];
    }
  }
  for (unsigned i = 0; i <= w; ++i) {
    std::vector<std::vector<int>> combs = combinations(k, i);
    for (const auto &comb : combs) {
      memcpy(c_temp, c_hard, n * sizeof(int));
      for (auto pos : comb) {
        for(unsigned int j=0; j<n; ++j) {
          c_temp[j] ^= g_buf[pos][j];
        }
      }
      double metr = 0;
      for (unsigned int i = 0; i < n; ++i) {
        metr += metric(c_temp[i], i);
        ++op_add;
      }
      if (min_metric != -1) {
        if (min_metric > metr) {
          min_metric = metr;
          memcpy(c_best, c_temp, n * sizeof(int));
        }
        ++op_cmp;
      } else {
        min_metric = metr;
        memcpy(c_best, c_temp, n * sizeof(int));
      }
    }
  }
  for (unsigned int i = 0; i < n; ++i) {
    u[permutation[i]] = c_best[i];
  }
  return min_metric;
}

void OrderedStatisticsDecoder::set_max_weight(int w) {
  this->w = w;
}

