#include <cmath>
#include <cstring>
#include "BSDDecoder.h"

void BSDDecoder::bubble_insert(Node *array, unsigned start, unsigned end, Node elem) {
  if (start == end) {
    array[start] = elem;
  }
  unsigned i;
  for (i = start; i < end; ++i) {
    if (elem.metric < array[i].metric) {
      for (unsigned j = i; j < end + 1; ++j) {
        std::swap(elem, array[j]);
      }
      break;
    }
    ++op_cmp;
  }
  if (i == end) {
    array[i] = elem;
  }
}

BSDDecoder::BSDDecoder(unsigned int n, unsigned int k, std::ifstream &filename) : TrellisDecoder(n, k, filename) {
  id = DecoderID::BSD;
  fwd_tree = new Node[2 * trellis_size];
  bkw_tree = new Node[2 * trellis_size];
}

BSDDecoder::BSDDecoder(unsigned int n, unsigned int k, int **h) : TrellisDecoder(n, k, h) {
  id = DecoderID::BSD;
  fwd_tree = new Node[2 * trellis_size];
  bkw_tree = new Node[2 * trellis_size];
}

BSDDecoder::~BSDDecoder() {
  delete[] fwd_tree;
  delete[] bkw_tree;
}

double BSDDecoder::decode(const double *x, int *u) {
  op_add = 0;
  op_cmp = 0;

  for (unsigned int i = 0; i < n; ++i) {
    alpha[i] = x[i] < 0 ? 0 : 1;
    ++op_cmp;
    beta[i] = fabs(x[i]);
  }
  memset(fwd_tree, 0, 2 * trellis_size * sizeof(Node));
  memset(bkw_tree, 0, 2 * trellis_size * sizeof(Node));

  double min_metric = -1;
  uint64_t min_candidate = 0;
  Node temp{};
  temp.tree = FWD;
  temp.layer = 0;
  temp.number = 0;
  temp.metric = 0;
  temp.path = 0;
  temp.path_avalaible[0] = true;
  temp.path_avalaible[1] = true;
  fwd_tree[0] = temp;
  temp.tree = BKW;
  temp.layer = n;
  bkw_tree[0] = temp;

  int fwdLayerReached = 0;
  int bkwLayerReached = n;
  int fwdTreeStart = 0;
  int fwdTreeEnd = 1;
  int bkwTreeStart = 0;
  int bkwTreeEnd = 1;

  uint64_t layerMask;
  while (min_metric == -1) {
    Node iter;
    if (fwdTreeEnd - fwdTreeStart > 0) {
      iter = fwd_tree[fwdTreeStart];
      fwd_tree[fwdTreeStart].tree = NIL;
      if (iter.layer < n) {
        layerMask = 0;
        for (unsigned int i = 0; i < k; ++i) {
          if (iter.layer < ranges[2 * i] || iter.layer >= ranges[2 * i + 1]) {
            layerMask ^= (uint64_t(1) << i);
          }
        }
        // Checking avalaible continuations for a node
        if (iter.number & layerMask) {
          iter.path_avalaible[0] = false;
        }
        uint64_t tempnum = iter.number;
        for (unsigned int i = 0; i < k; ++i) {
          tempnum ^= (uint64_t(h[i][iter.layer]) << i);
        }
        if (tempnum & layerMask) {
          iter.path_avalaible[1] = false;
        }
        for (unsigned int j = 0; j < 2; ++j) {
          // If path is not available, pathAvailable <- false and skipping
          if (iter.path_avalaible[j]) {
            temp.tree = FWD;
            temp.layer = iter.layer + 1;
            fwdLayerReached = fwdLayerReached < temp.layer ? temp.layer : fwdLayerReached;
            temp.number = j ? tempnum : iter.number;
            double m = metric(j, iter.layer);
            temp.path = j ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
            temp.path_avalaible[0] = true;
            temp.path_avalaible[1] = true;
            if (m != 0) {
              temp.metric = iter.metric + m;
              ++op_add;
              bubble_insert(fwd_tree, fwdTreeStart + 1, fwdTreeEnd, temp);
              ++fwdTreeEnd;
            } else {
              temp.metric = iter.metric;
              fwd_tree[fwdTreeStart] = temp;
            }
            ++op_cmp;
          }
        }
        if (fwd_tree[fwdTreeStart].tree == NIL && fwdTreeStart < fwdTreeEnd) {
          ++fwdTreeStart;
        }
      }
    }
    if (bkwLayerReached - fwdLayerReached <= 0) {
      for (unsigned l = bkwTreeStart; l < bkwTreeEnd; ++l) {
        for (unsigned j = fwdTreeStart; j < fwdTreeEnd; ++j) {
          if (fwd_tree[j].layer == bkw_tree[l].layer &&
              fwd_tree[j].number == bkw_tree[l].number) {
            double m = fwd_tree[j].metric + bkw_tree[l].metric;
            ++op_add;
            if (min_metric == -1) {
              min_metric = m;
              min_candidate = bkw_tree[l].path + fwd_tree[j].path;
            } else {
              if (m < min_metric) {
                min_metric = m;
                min_candidate = bkw_tree[l].path + fwd_tree[j].path;
              }
              ++op_cmp;
            }
          }
        }
      }
    }
    if (bkwTreeEnd - bkwTreeStart > 0) {
      iter = bkw_tree[bkwTreeStart];
      bkw_tree[bkwTreeStart].tree = NIL;
      if (iter.layer > 0) {
        layerMask = 0;
        for (unsigned int i = 0; i < k; ++i) {
          if ((iter.layer - 1) <= ranges[2 * i] || (iter.layer - 1) > ranges[2 * i + 1]) {
            layerMask ^= (uint64_t(1) << i);
          }
        }

        if (iter.number & layerMask) {
          iter.path_avalaible[0] = false;
        }
        uint64_t tempnum = iter.number;
        for (unsigned int i = 0; i < k; ++i) {
          tempnum ^= (uint64_t(h[i][iter.layer - 1]) << i);
        }
        if (tempnum & layerMask) {
          iter.path_avalaible[1] = false;
        }
        for (unsigned int j = 0; j < 2; ++j) {
          // If path is not available, pathAvailable <- false and skipping
          if (iter.path_avalaible[j]) {
            temp.tree = BKW;
            temp.layer = iter.layer - 1;
            bkwLayerReached = bkwLayerReached > temp.layer ? temp.layer : bkwLayerReached;
            temp.number = j ? tempnum : iter.number;
            double m = metric(j, iter.layer - 1);
            temp.path = j ? iter.path ^ (uint64_t(1) << (iter.layer - 1)) : iter.path;
            temp.path_avalaible[0] = true;
            temp.path_avalaible[1] = true;
            if (m != 0) {
              temp.metric = iter.metric + m;
              ++op_add;
              bubble_insert(bkw_tree, bkwTreeStart + 1, bkwTreeEnd, temp);
              ++bkwTreeEnd;
            } else {
              temp.metric = iter.metric;
              bkw_tree[bkwTreeStart] = temp;
            }
            ++op_cmp;
          }
        }
        if (bkw_tree[bkwTreeStart].tree == NIL && bkwTreeStart < bkwTreeEnd) {
          ++bkwTreeStart;
        }
      }
    }
    if (bkwLayerReached - fwdLayerReached <= 0) {
      for (unsigned l = fwdTreeStart; l < fwdTreeEnd; ++l) {
        for (unsigned j = bkwTreeStart; j < bkwTreeEnd; ++j) {
          if (bkw_tree[j].layer == fwd_tree[l].layer &&
              fwd_tree[l].number == bkw_tree[j].number) {
            double m = fwd_tree[l].metric + bkw_tree[j].metric;
            ++op_add;
            if (min_metric == -1) {
              min_metric = m;
              min_candidate = bkw_tree[j].path + fwd_tree[l].path;
            } else {
              if (m < min_metric) {
                min_metric = m;
                min_candidate = bkw_tree[j].path + fwd_tree[l].path;
              }
              ++op_cmp;
            }
          }
        }
      }
    }
  }
  for (unsigned int i = 0; i < n; ++i) {
    u[i] = (min_candidate & (uint64_t(1) << i)) ? 1 : 0;
  }
  return min_metric;
}
