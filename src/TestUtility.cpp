#include "TestUtility.h"

void exhaustive_subtrellis_verification(unsigned int n, unsigned int k, int** g, Trellis &trellis, unsigned int w, int** check) {
  int *codeword = new int[n];
  unsigned int word_weight = 0;
  for(uint64_t infoword=1; infoword<(uint64_t(1)<< k); ++infoword) {
    memset(codeword, 0, n*sizeof(int));
    word_weight = 0;
    for(unsigned int j=0; j<n; ++j) {
      int temp = 0;
      for(unsigned int l=0; l<k; ++l) {
        temp ^= g[l][j] & ((infoword & (uint64_t(1) << l)) ? 1 : 0);
      }
      word_weight += temp ? 1 : 0;
      codeword[j] = temp ? 1 : 0;
    }
    if(check != nullptr) {
      unsigned int synd_weight = 0;
      for (unsigned int j = 0; j < n - k; ++j) {
        int temp = 0;
        for (unsigned int i = 0; i < n; ++i) {
          temp ^= codeword[i] & check[j][i];
        }
        synd_weight += temp;
      }
      INFO("Infoword: " << std::bitset<64>(infoword).to_string());
      INFO("Codeword: " << arrayToSstream<int>(n, codeword).str());
      REQUIRE(synd_weight == 0);
    }
    bool should_be_present = word_weight <= w;
    bool present = true;
    unsigned long long state = 0;
    for(unsigned int j=0; j<n; ++j) {
      state = trellis[j][state].next_node[codeword[j]];
      if(state == ~0) {
        present = false;
        break;
      }
    }
    INFO("Infoword: " << std::bitset<64>(infoword).to_string());
    INFO("Codeword: " << arrayToSstream<int>(n, codeword).str());
    INFO("Weight: " << word_weight);
    INFO("Present: " << present);
    INFO("Should be present: " <<should_be_present);
    CHECK(((present && should_be_present) || (!present && !should_be_present)));
  }
  delete[] codeword;
}

void test_beast_decoder(int tests, double stn, double delta, const std::string &filename) {
  std::ifstream in(filename);
  Simulation sim(in, 0, tests);
  sim.add_beast_decoder(delta);
  sim.setSTN(stn);
  sim.test_run();
}
void test_bsd_decoder(int tests, double stn, std::string filename) {
  std::ifstream in(filename);
  Simulation sim(in, 0, tests);
  sim.add_bsd_decoder();
  sim.setSTN(stn);
  sim.test_run();
}

void test_KTKL_decoder(int tests, double stn, int w, int buf_size, const std::string &filename) {
  std::ifstream in(filename);
  Simulation sim(in, 0, tests);
  sim.add_ktkl_decoder(w, buf_size);
  sim.setSTN(stn);
  sim.test_run();
}
void test_osd_decoder(int tests, double stn, int w, std::string filename) {
  std::ifstream in(filename);
  Simulation sim(in, 0, tests);
  sim.add_ordered_statistics_decoder(w);
  sim.setSTN(stn);
  sim.test_run();
}
void test_viterbi_decoder(int tests, double stn, const std::string &filename) {
  std::ifstream in(filename);
  Simulation sim(in, 0, tests);
  sim.add_viterbi_decoder();
  sim.setSTN(stn);
  sim.test_run();
}
