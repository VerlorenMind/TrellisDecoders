
#include <Trellis.h>
#include <Utility.h>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <iostream>
#include <sstream>

TrellisNode &TrellisLayer::operator[](unsigned int i) {
  return layer[i];
}

void TrellisLayer::init(unsigned int size, unsigned int num) {
  layer_size = size;
  layer_num = num;
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

TrellisLayer &Trellis::operator[](unsigned int i) {
  return trellis[i];
}

Trellis::Trellis() {
  trellis = nullptr;
  max_layer_size = 0;
}

/*
void Trellis::construct_from_check_matrix(unsigned int n, unsigned int k, int **h) {
  trellis = new TrellisLayer[n+1];
  trellis[0].init(1, 0);
  minspanForm(n, k, h);
  unsigned int *ranges = findRanges(n, k , h);
  unsigned int *active_bits = new unsigned int[k];
  unsigned int num_of_act_bits = 0;
  for(unsigned int i=1; i<=n; ++i) {
    for(unsigned int j=0; j<k; ++j) {
      if(ranges[2*j] < i && ranges[2*j+1]) {
        active_bits[num_of_act_bits++] = j;
      }
    }
    trellis[i].init((uint64_t(1) << num_of_act_bits), i);
    for(unsigned int j=0; j<trellis[i-1].size(); ++j) {

    }
  }
}
*/
void Trellis::construct_from_gen_matrix(unsigned int n, unsigned int k, int **gen) {
  trellis = new TrellisLayer[n + 1];
  std::vector<unsigned> row_start(n);
  std::vector<unsigned> row_end(n);
  unsigned *active_bits = new unsigned[k]; //the list of active input bits
  uint64_t *cw_0 = new uint64_t[1ull << (n / 2 + 1)];
  uint64_t *cw_1 = new uint64_t[1ull << (n / 2 + 1)];
  unsigned *num_of_active_bits = new unsigned[n + 1];
  int **g = new int *[k];
  for (unsigned int i = 0; i < k; ++i) {
    g[i] = new int[n];
    memcpy(g[i], gen[i], n * sizeof(int));
  }
  minspanForm(n, k, g, row_start, row_end);
  std::string tempstr = matrixToSstream(k, n, g).str();
  trellis = new TrellisLayer[n + 1];
  unsigned num_of_AB = 0, next_num_of_AB = 0;
  cw_0[0] = 0;
  num_of_active_bits[0] = 0;
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
    if (j < n) {
      num_of_active_bits[j + 1] = num_of_AB;
    }
    if (num_of_AB > max_layer_size)
      max_layer_size = num_of_AB;
  };

  std::stringstream temp;
  temp << "digraph G{\n";
  for (unsigned j = 0; j < n; j++) {
    for (unsigned s = 0; s < 1ull << num_of_active_bits[j]; s++) {
      if (trellis[j][s].next_node[0] != ~0)
        temp << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_'
                  << trellis[j][s].next_node[0] << "[style=dashed];\n";
      if (trellis[j][s].next_node[1] != ~0)
        temp << "S_" << j << '_' << s << "->" << "S_" << (j + 1) << '_'
                  << trellis[j][s].next_node[1] << "[style=solid];\n";
    }
  }
  temp << "}\n===========\n";
  tempstr = temp.str();

  max_layer_size = 1ull << max_layer_size;
  delete[] cw_0;
  delete[] cw_1;
  delete[] active_bits;
  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
  delete[] num_of_active_bits;
}
unsigned long long Trellis::get_max_layer_size() {
  return max_layer_size;
}
Trellis::~Trellis() {
  delete[] trellis;
}

