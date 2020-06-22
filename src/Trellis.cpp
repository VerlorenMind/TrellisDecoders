
#include <Trellis.h>
#include <Utility.h>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <iostream>
#include <sstream>
#include <stack>
#include <unordered_set>

TrellisNode &TrellisLayer::operator[](unsigned int i) {
  return layer[i];
}

void TrellisLayer::init(unsigned int size, unsigned int num) {
  layer_size = size;
  true_size = size;
  layer = new TrellisNode[size];
  for (unsigned int i = 0; i < size; ++i) {
    layer[i].next_node[0] = ~0;
    layer[i].next_node[1] = ~0;
    layer[i].prev_node[0] = ~0;
    layer[i].prev_node[1] = ~0;
  }
}
unsigned TrellisLayer::size() {
  return layer_size;
}
TrellisLayer::~TrellisLayer() {
  delete[] layer;
}
void TrellisLayer::delete_node(unsigned long long num) {
  memmove(layer + num, layer + num + 1, (layer_size - num - 1) * sizeof(TrellisNode));
  --layer_size;
}
void TrellisLayer::shrink_to_fit() {
  if (true_size > layer_size) {
    TrellisNode *new_layer = new TrellisNode[layer_size];
    memcpy(new_layer, layer, layer_size * sizeof(TrellisNode));
    delete[] layer;
    layer = new_layer;
    true_size = layer_size;
  }
}
void TrellisLayer::add_node() {
  if(layer_size == true_size) {
    TrellisNode *new_layer = new TrellisNode[layer_size + 1];
    true_size = layer_size + 1;
    memcpy(new_layer, layer, layer_size*sizeof(TrellisNode));
    delete[] layer;
    layer = new_layer;
  }
  for(unsigned z=0; z<2; ++z) {
    layer[layer_size].next_node[z] = ~0;
    layer[layer_size].prev_node[z] = ~0;
  }
  ++layer_size;
}

TrellisLayer &Trellis::operator[](unsigned int i) {
  return trellis[i];
}

Trellis::Trellis() {
  trellis = nullptr;
  max_layer_size = 0;
}

void Trellis::construct_from_gen_matrix(unsigned int n, unsigned int k, int **gen) {
  trellis_size = n + 1;
  trellis = new TrellisLayer[n + 1];
  std::vector<unsigned> row_start(n);
  std::vector<unsigned> row_end(n);
  unsigned *active_bits = new unsigned[n]; //the list of active input bits
  uint64_t *cw_0 = new uint64_t[1ull << (n / 2 + 1)];
  uint64_t *cw_1 = new uint64_t[1ull << (n / 2 + 1)];
  // unsigned *num_of_active_bits = new unsigned[n + 1];
  int **g = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    g[i] = new int[n];
    memcpy(g[i], gen[i], n * sizeof(int));
  }
  minspanForm(n, k, g, row_start, row_end);
  std::string tempstr = matrixToSstream(k, n, g).str();
  unsigned num_of_AB = 0, next_num_of_AB = 0;
  cw_0[0] = 0;
  // num_of_active_bits[0] = 0;
  trellis[0].init(1ull << num_of_AB, 0);
  for (unsigned j = 0; j < n; j++) {
    size_t b = std::find(active_bits, active_bits + num_of_AB, row_end[j]) - active_bits;
    unsigned long long extraction_mask = (row_end[j] == ~0) ? ~0ull : ((1ull << b) - 1);
    if (row_start[j] != ~0) {
      active_bits[next_num_of_AB++] = row_start[j];
    }
    if (row_end[j] != ~0) {
      memmove(active_bits + b, active_bits + b + 1, sizeof(unsigned) * (next_num_of_AB - b));
      --next_num_of_AB;
    };

    if (j < n) {
      trellis[j + 1].init(1ull << next_num_of_AB, j + 1);
      for (unsigned long long state = 0; state < 1ull << next_num_of_AB; state++) {
        for (int k = 0; k < 2; ++k) {
          trellis[j + 1][state].prev_node[k] = (unsigned long long int) ~0;
        }
      }
    }
    if (row_start[j] == ~0) {
      //no new rows are starting here
      for (unsigned long long state = 0; state < 1ull << num_of_AB; state++) {
        uint64_t next_bit = (cw_0[state] >> j) & 1;
        unsigned long long next_state = (state & extraction_mask) | ((state >> 1) & ~extraction_mask);
        cw_1[next_state] = cw_0[state];
        trellis[j][state].next_node[next_bit] = next_state;
        trellis[j][state].next_node[1 - next_bit] = ~0;
        if (j < n) {
          trellis[j + 1][next_state].prev_node[next_bit] = state;
        }
      };
    } else {

      for (unsigned long long state = 0; state < 1ull << num_of_AB; state++) {
        unsigned long long next_state_0 = state;
        unsigned long long next_state_1 = state ^(1ull << num_of_AB);
        //compress the state variables
        next_state_0 = (next_state_0 & extraction_mask) | ((next_state_0 >> 1) & ~extraction_mask);
        next_state_1 = (next_state_1 & extraction_mask) | ((next_state_1 >> 1) & ~extraction_mask);
        uint64_t c1 = cw_0[state];
        for (unsigned int i = 0; i < n; ++i) {
          c1 ^= (g[row_start[j]][i] << i);
        }
        cw_1[next_state_0] = cw_0[state];
        cw_1[next_state_1] = c1;
        uint64_t next_bit_0 = (cw_0[state] >> j) & 1;
        uint64_t next_bit_1 = (c1 >> j) & 1;
        assert((next_bit_0 ^ next_bit_1) == 1);
        trellis[j][state].next_node[next_bit_0] = next_state_0;
        trellis[j][state].next_node[next_bit_1] = next_state_1;
        if (j < n) {
          trellis[j + 1][next_state_0].prev_node[next_bit_0] = state;
          trellis[j + 1][next_state_1].prev_node[next_bit_1] = state;
        }
      }
    };
    std::swap(cw_0, cw_1);
    num_of_AB = next_num_of_AB;
    /*if (j < n) {
      num_of_active_bits[j + 1] = num_of_AB;
    }*/
    if (num_of_AB > max_layer_size)
      max_layer_size = num_of_AB;
  };

  max_layer_size = 1ull << max_layer_size;
  delete[] cw_0;
  delete[] cw_1;
  delete[] active_bits;
  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
  // delete[] num_of_active_bits;
}
unsigned long long Trellis::get_max_layer_size() {
  return max_layer_size;
}
Trellis::~Trellis() {
  delete[] trellis;
}
void Trellis::print_trellis(std::ostream &out) {
  out << "digraph G{\n";
  for (unsigned j = 0; j < trellis_size; j++) {
    for (unsigned s = 0; s < trellis[j].size(); s++) {
      if (trellis[j][s].next_node[0] != ~0)
        out << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_'
            << trellis[j][s].next_node[0] << "[style=dashed];\n";
      if (trellis[j][s].next_node[1] != ~0)
        out << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_'
            << trellis[j][s].next_node[1] << "[style=solid];\n";
    }
  }
  out << "}";
}
void Trellis::delete_node(unsigned int layer, unsigned long long num) {
  for (unsigned int z = 0; z < 2; ++z) {
    if (trellis[layer][num].prev_node[z] != ~0) {
      trellis[layer - 1][trellis[layer][num].prev_node[z]].next_node[z] = ~0;
    }
    if (trellis[layer][num].next_node[z] != ~0) {
      trellis[layer + 1][trellis[layer][num].next_node[z]].prev_node[z] = ~0;
    }
  }
  trellis[layer].delete_node(num);
  for (unsigned long t = num; t < trellis[layer].size(); ++t) {
    for (unsigned int z = 0; z < 2; ++z) {
      if (trellis[layer][t].prev_node[z] != ~0) {
        trellis[layer - 1][trellis[layer][t].prev_node[z]].next_node[z] = t;
      }
      if (trellis[layer][t].next_node[z] != ~0) {
        trellis[layer + 1][trellis[layer][t].next_node[z]].prev_node[z] = t;
      }
    }
  }
  if(layer == trellis_size - 2 && num == second_0_edge_to_end) {
    second_0_edge_to_end = ~0;
  }
}
bool Trellis::dfs_trellis_search(std::pair<unsigned int, unsigned long long> node,
                                 std::vector<bool> *keep_node, unsigned int metric, unsigned int w,
                                 std::vector<int> &path, std::vector<std::vector<int>> &paths_to_keep) {
  if (metric > w) {
    return false;
  }
  if (node.first == (trellis_size - 1)) {
    keep_node[node.first][node.second] = true;
    paths_to_keep.push_back(path);
    return true;
  }
  bool child[2];
  child[0] = false;
  child[1] = false;
  for (unsigned int z = 0; z < 2; ++z) {
    if (trellis[node.first][node.second].next_node[z] != ~0) {
      path[node.first] = z;
      child[z] = dfs_trellis_search(std::make_pair(node.first + 1, trellis[node.first][node.second].next_node[z]),
                                    keep_node,
                                    metric + z,
                                    w, path, paths_to_keep);
    }
  }
  if (child[0] || child[1]) {
    keep_node[node.first][node.second] = true;
    return true;
  } else {
    return false;
  }
}

void Trellis::reduce_to_weight_dfs(unsigned int w) {
  std::vector<bool> *keep_node = new std::vector<bool>[trellis_size];
  std::vector<std::vector<int>> words_to_keep;
  std::vector<int> word(trellis_size - 1);
  for (unsigned int i = 0; i < trellis_size; ++i) {
    keep_node[i].resize(trellis[i].size());
    std::fill(keep_node[i].begin(), keep_node[i].end(), false);
  }
  dfs_trellis_search(std::make_pair(0, 0), keep_node, 0, w, word, words_to_keep);
  TrellisLayer *new_trellis = new TrellisLayer[trellis_size];
  for (unsigned int i = 0; i < trellis_size; ++i) {
    unsigned long long trellis_state = 0;
    for (unsigned long long j = 0; j < trellis[i].size(); ++j) {
      if(!keep_node[i][j]) {
        delete_node(i, trellis_state);
      }
      else {
        ++trellis_state;
      }
    }
    trellis[i].shrink_to_fit();
    new_trellis[i].init(trellis[i].size(), i);
  }
  unsigned long long state;
/*
  std::ofstream out;
  out.open("../tests/trellis.gv");
  print_trellis(out);
  out.close();
  system("dot ../tests/trellis.gv -Tpng -o ../tests/inter2.png");
*/
  for (auto kept_word : words_to_keep) {
    state = 0;
    for (unsigned int i = 0; i < trellis_size - 1; ++i) {
      new_trellis[i][state].next_node[kept_word[i]] = trellis[i][state].next_node[kept_word[i]];
      state = new_trellis[i][state].next_node[kept_word[i]];
    }
  }
  delete[] trellis;
  trellis = new_trellis;
  delete[] keep_node;
}

void Trellis::reduce_to_weight(unsigned int w) {
  unsigned int **min_weight_to, **min_weight_from;
  min_weight_to = new unsigned int*[trellis_size];
  min_weight_from = new unsigned int*[trellis_size];
  for(unsigned int i=0; i<trellis_size; ++i) {
    min_weight_to[i] = new unsigned int[trellis[i].size()];
    min_weight_from[i] = new unsigned int[trellis[i].size()];
    // Filling the arrays with maximum values
    memset(min_weight_to[i], ~0, trellis[i].size() * sizeof(unsigned int));
    memset(min_weight_from[i], ~0, trellis[i].size() * sizeof(unsigned int));
  }
  min_weight_to[0][0] = 0;
  min_weight_from[trellis_size-1][0] = 0;
  // Finding minimum path weight to all nodes from the first
  for (unsigned j = 0; j < trellis_size - 1; j++) {
    for (unsigned l = 0; l < trellis[j].size(); l++) {
      for (unsigned z = 0; z < 2; z++) {
        unsigned long long s1 = trellis[j][l].next_node[z];
        if (s1 != ~0) {
          unsigned int metric = min_weight_to[j][l] + z;
          if (metric < min_weight_to[j+1][s1]) {
            min_weight_to[j+1][s1] = metric;
          }
        }
      }
    }
  }
  // Finding minimum path weight from all nodes to the final
  for (unsigned j = trellis_size - 1; j > 0; j--) {
    for (unsigned l = 0; l < trellis[j].size(); l++) {
      for (unsigned z = 0; z < 2; z++) {
        unsigned long long s1 = trellis[j][l].prev_node[z];
        if (s1 != ~0) {
          unsigned int metric = min_weight_from[j][l] + z;
          if (metric < min_weight_from[j-1][s1]) {
            min_weight_from[j-1][s1] = metric;
          }
        }
      }
    }
  }
  // Purging trellis
  for(unsigned i = 1; i < trellis_size; ++i) {
    unsigned long long state = 0;
    unsigned long long original_size = trellis[i].size();
    for(unsigned long long l = 0; l < original_size; ++l) {
      // Purging if the weight of the minimum path is greater than w
      if(min_weight_to[i][l] + min_weight_from[i][l] > w) {
        delete_node(i, state);
      }
      // If no incoming edges, purge
      else if(trellis[i][state].prev_node[0] == ~0 && trellis[i][state].prev_node[1] == ~0) {
        delete_node(i, state);
      }
      // Purging branches of paths with excessive weight
      else {
        unsigned long long next_state;
        for(unsigned z=0; z<2; ++z) {
          next_state = trellis[i][state].next_node[z];
          if(next_state != ~0 && min_weight_to[i][l] + min_weight_from[i+1][next_state] + z > w) {
            trellis[i][state].next_node[z] = ~0;
            trellis[i+1][next_state].prev_node[z] = ~0;
          }
        }
        ++state;
      }
    }
  }
  // Dealing with leftover paths of weight above w
  // Finding the first non-zero merging into all-zero path
  unsigned start_index = 1;
  for(; start_index<trellis_size-1; ++start_index) {
    if(trellis[start_index][0].prev_node[1] != ~0) {
      trellis[start_index].shrink_to_fit();
      break;
    }
  }
  // Adding additional nodes to reroute the mergings
  if(start_index < trellis_size - 1) {
    for (unsigned i = start_index; i < trellis_size - 1; ++i) {
      trellis[i].add_node();
      if (i > start_index) {
        trellis[i][trellis[i].size() - 1].prev_node[0] = trellis[i - 1].size() - 1;
        trellis[i - 1][trellis[i - 1].size() - 1].next_node[0] = trellis[i].size() - 1;
      }
      unsigned long long merging_state = trellis[i][0].prev_node[1];
      if (merging_state != ~0) {
        trellis[i][0].prev_node[1] = ~0;
        trellis[i][trellis[i].size() - 1].prev_node[1] = merging_state;
        trellis[i - 1][merging_state].next_node[1] = trellis[i].size() - 1;
      }
    }
    trellis[trellis_size - 2][trellis[trellis_size - 2].size() - 1].next_node[0] = 0;
    second_0_edge_to_end = trellis[trellis_size - 2].size() - 1;
  }
  // Deleting the all-zero word
  for(unsigned i = trellis_size-2; i >= 0 && trellis[i][0].next_node[1] == ~0; --i) {
    delete_node(i, 0);
    trellis[i].shrink_to_fit();
  }
  for(unsigned i = 0; i<trellis_size; ++i) {
    delete[] min_weight_to[i];
    delete[] min_weight_from[i];
  }
  delete[] min_weight_to;
  delete[] min_weight_from;
}
void Trellis::construct_from_check_matrix(unsigned int n, unsigned int k, int **check) {
  trellis_size = n + 1;
  trellis = new TrellisLayer[n + 1];
  int **h = new int *[n-k];
  for (unsigned int i = 0; i < n-k; ++i) {
    h[i] = new int[n];
    memcpy(h[i],check[i], n * sizeof(int));
  }
  std::vector<unsigned> row_start(n);
  std::vector<unsigned> row_end(n);
  minspanForm(n, n-k, h, row_start, row_end);
  std::vector<bool> active_bits_mask(n-k);
  unsigned int active_bits = 0;
  std::vector<std::vector<bool>> prev_accumulated_synd(1), cur_accumulated_synd;
  trellis[0].init(1, 0);
  prev_accumulated_synd[0].resize(n-k);
  std::fill(prev_accumulated_synd[0].begin(), prev_accumulated_synd[0].end(), 0);
  for(unsigned int i=0; i<n; ++i) {
    if(row_start[i] != ~0) {
      active_bits_mask[row_start[i]] = true;
      ++active_bits;
    }
    if(row_end[i] != ~0) {
      --active_bits;
      active_bits_mask[row_end[i]] = false;
    }
    trellis[i+1].init(uint64_t(1) << active_bits, i+1);
    for(unsigned int j=0; j<trellis[i].size(); ++j) {
      std::vector<bool> synd(n-k);
      std::copy(prev_accumulated_synd[j].begin(), prev_accumulated_synd[j].end(), synd.begin());
      for(unsigned int z=0; z<2; ++z) {
        if(z) {
          for(unsigned int l=0; l<n-k; ++l) {
            synd[l] = (h[l][i] && !synd[l]) || (!h[l][i] && synd[l]);
          }
        }
        bool valid = true;
        for (unsigned int l = 0; l < n - k; ++l) {
          if (synd[l] && !active_bits_mask[l]) {
            valid = false;
            break;
          }
        }
        if (valid) {
          bool found = false;
          for (unsigned long long l = 0; l < cur_accumulated_synd.size(); ++l) {
            std::vector<bool> &accum_synd = cur_accumulated_synd[l];
            bool equal = true;
            for(unsigned int b=0; b<n-k; ++b) {
              if(accum_synd[b] != synd[b]) {
                equal = false;
                break;
              }
            }
            if (equal) {
              found = true;
              trellis[i][j].next_node[z] = l;
              trellis[i+1][l].prev_node[z] = j;
              break;
            }
          }
          if (!found) {
            cur_accumulated_synd.push_back(synd);
            trellis[i][j].next_node[z] = cur_accumulated_synd.size() - 1;
            trellis[i+1][cur_accumulated_synd.size()-1].prev_node[z] = j;
          }
        }
      }
    }
    prev_accumulated_synd.swap(cur_accumulated_synd);
    cur_accumulated_synd.clear();
  }
  for(unsigned int i=0; i<n-k; ++i) {
    delete[] h[i];
  }
  delete[] h;
}

