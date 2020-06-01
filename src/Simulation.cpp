#include <fstream>
#include <cmath>
#include <random>
#include <Utility.h>
#include <BeastDecoder.h>
#include <BSDDecoder.h>
#include <iostream>
#include <chrono>
#include <OrderedStatisticsDecoder.h>
#include <ViterbiDecoder.h>
#include "Simulation.h"

const double EPSILON = 0.0000001;

void Simulation::generate_information_vector() {
  static std::uniform_int_distribution<int> distribution(0, 1);
  for (unsigned int i = 0; i < k; ++i) {
    u[i] = distribution(gen);
  }
}

void Simulation::encode() {
  int temp = 0;
  for (unsigned int i = 0; i < n; ++i) {
    temp = 0;
    for (unsigned int j = 0; j < k; ++j) {
      temp ^= g[j][i] & u[j];
    }
    ux[i] = (temp ? 1 : 0);
  }
}

void Simulation::apply_noise() {
  // Applying AWGN with BPSK modulation
  std::normal_distribution<double> rv(0.0, dev);
  for (unsigned int i = 0; i < n; ++i) {
    x[i] = (ux[i] ? 1.0 : -1.0) + rv(gen);
  }
}

int Simulation::calculate_syndrome(const int *vec) {
  memset(synd, 0, (n - k) * sizeof(int));
  int syndsum = 0;
  for (unsigned int j = 0; j < n - k; ++j) {
    for (unsigned int i = 0; i < n; ++i) {
      synd[j] ^= vec[i] & h[j][i];
    }
  }
  for (unsigned int j = 0; j < n - k; ++j) {
    syndsum += synd[j];
  }
  return syndsum;
}

bool Simulation::check_decoded_word() {
  bool flag = true;
  for (unsigned int i = 0; i < n; ++i) {
    if (y[i] != ux[i]) {
      flag = false;
      break;
    }
  }
  return flag;
}

Simulation::Simulation(std::ifstream &input_file, unsigned int max_errors,
                       unsigned int max_iterations, unsigned int seed) {
  std::getline(input_file, code_name);
  // input_file >> code_name;
  input_file >> n >> k;
  g = readMatrix(input_file, n, k);
  h = readMatrix(input_file, n, n - k);
  test = 0;
  this->max_iterations = max_iterations;
  this->max_errors = max_errors;
  dev = 0;
  stn = 0;
  if (seed == 0) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  gen.seed(seed);
  x = new double[n];
  y = new int[n];
  u = new int[k];
  ux = new int[n];
  synd = new int[n - k];
}

Simulation::Simulation(std::ifstream &input_file, std::ifstream &run_config) : Simulation(input_file) {
  setSTN(1);
  std::string line;
  std::getline(run_config, line);
  std::stringstream linestr(line);
  linestr >> max_errors >> max_iterations;
  unsigned int seed;
  linestr >> seed;
  if (seed == 0) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  gen.seed(seed);
  decoders.resize(0);
  std::string decoder_name;
  while (std::getline(run_config, line)) {
    linestr.str(line);
    linestr.clear();
    linestr >> decoder_name;
    DecoderID id = stringToId(decoder_name);
    switch (id) {
      case DecoderID::BEAST: {
        double delta;
        linestr >> delta;
        SoftDecoder *beast = new BeastDecoder(n, k, h, delta);
        decoders.push_back(beast);
        break;
      }
      case DecoderID::BSD: {
        SoftDecoder *bsd = new BSDDecoder(n, k, h);
        decoders.push_back(bsd);
        break;
      }
      case DecoderID::KTKL: {
        break;
      }
      case DecoderID::ERROR: {
        std::cerr << "Invalid run config: " + decoder_name + " is not a name";
      }
      case DecoderID::ORDERED_STATISTICS: {
        int w;
        linestr >> w;
        SoftDecoder *osd = new OrderedStatisticsDecoder(n, k, g, w);
        decoders.push_back(osd);
      }
      case DecoderID::VITERBI: {
        SoftDecoder *vit = new ViterbiDecoder(n, k, g);
        decoders.push_back(vit);
      }
    }
  }
  errors_by_decoder.resize(decoders.size());
  cmp_ops_by_decoder.resize(decoders.size());
  add_ops_by_decoder.resize(decoders.size());
}

Simulation::Simulation(std::ifstream &input_file, double stn, SoftDecoder *decoder, unsigned int max_errors,
                       unsigned int max_iterations, unsigned int seed) : Simulation(input_file) {
  setSTN(stn);
  if (seed == 0) {
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  }
  gen.seed(seed);
  this->max_errors = max_errors;
  this->max_iterations = max_iterations;
  decoders.resize(1);
  decoders[0] = decoder;
  errors_by_decoder.resize(decoders.size());
  cmp_ops_by_decoder.resize(decoders.size());
  add_ops_by_decoder.resize(decoders.size());
}

Simulation::~Simulation() {
  delete[] x;
  delete[] y;
  delete[] u;
  for (unsigned int i = 0; i < k; ++i) {
    delete[] g[i];
  }
  delete[] g;
  for (unsigned int i = 0; i < n - k; ++i) {
    delete[] h[i];
  }
  delete[] h;
  delete[] ux;
  delete[] synd;
  for (unsigned int i = 0; i < decoders.size(); ++i) {
    delete decoders[i];
  }
}
void Simulation::run() {
  errors_by_decoder.resize(decoders.size());
  std::fill(errors_by_decoder.begin(), errors_by_decoder.end(), 0);
  cmp_ops_by_decoder.resize(decoders.size());
  std::fill(cmp_ops_by_decoder.begin(), cmp_ops_by_decoder.end(), 0);
  add_ops_by_decoder.resize(decoders.size());
  std::fill(add_ops_by_decoder.begin(), add_ops_by_decoder.end(), 0);
  std::vector<bool> decoder_stopped(decoders.size());
  bool all_stopped = false;
  std::fill(decoder_stopped.begin(), decoder_stopped.end(), false);
  test = 0;
  while (!all_stopped && test < max_iterations) {
    generate_information_vector();
    encode();
    apply_noise();
    for (unsigned long i = 0; i < decoders.size(); ++i) {
      if (decoder_stopped[i]) {
        continue;
      }
      decoders[i]->decode(x, y);
      cmp_ops_by_decoder[i] += decoders[i]->op_cmp;
      add_ops_by_decoder[i] += decoders[i]->op_add;
      if (!check_decoded_word()) {
        ++errors_by_decoder[i];
        if (errors_by_decoder[i] >= max_errors) {
          decoder_stopped[i] = true;
        }
      }
    }
    all_stopped = true;
    for (bool stop : decoder_stopped) {
      all_stopped = all_stopped && stop;
    }
    ++test;
  }
}

SoftDecoder *Simulation::get_decoder(DecoderID id) {
  for (auto dec : decoders) {
    if (dec->get_id() == id) {
      return dec;
    }
  }
  return nullptr;
}

void Simulation::add_decoder(DecoderID id) {
  switch (id) {
    case DecoderID::BEAST: {
      SoftDecoder *beast = new BeastDecoder(n, k, h, 0.5);
      decoders.push_back(beast);
      break;
    }
    case DecoderID::BSD: {
      SoftDecoder *bsd = new BSDDecoder(n, k, h);
      decoders.push_back(bsd);
      break;
    }
    case DecoderID::ORDERED_STATISTICS: {
      SoftDecoder *osd = new OrderedStatisticsDecoder(n, k, g, 0);
      decoders.push_back(osd);
      break;
    }
    case DecoderID::VITERBI: {
      SoftDecoder *vit = new ViterbiDecoder(n, k, g);
      decoders.push_back(vit);
      break;
    }
    case DecoderID::KTKL: {
      return;
    }
    case DecoderID::ERROR: {
      return;
    }
  }
  errors_by_decoder.resize(decoders.size());
  cmp_ops_by_decoder.resize(decoders.size());
  add_ops_by_decoder.resize(decoders.size());
}
void Simulation::setSTN(double stn) {
  this->stn = stn;
  dev = sqrt((1 / (((1.0 * k) / n) * pow(10, stn / 10))) / 2);
}

#ifdef CATCH_TESTING
void Simulation::test_run() {
  std::normal_distribution<double> rv(0.0, dev);
  double metric, calcWeight, euclTrue, euclCalc, overall = 0;
  bool flag;
  std::string word;
  int fails = 0;
  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
  for (test = 0; test < max_iterations; ++test) {
    INFO("Test #" << test);
    generate_information_vector();
    INFO("Informational word: " << arrayToSstream<int>(k, u).str());
    encode();
    INFO("Coded word: " << arrayToSstream<int>(n, ux).str());
    apply_noise();
    INFO("Coded word syndrome: " << arrayToSstream<int>(n - k, synd).str());
    REQUIRE(calculate_syndrome(ux) == 0);
    INFO("Coded word with noise: " << arrayToSstream<double>(n, x).str());
    metric = 0;
    for (unsigned int i = 0; i < n; ++i) {
      if (x[i] < 0) {
        metric += (ux[i] == 0 ? 0 : fabs(x[i]));
      } else {
        metric += (ux[i] == 1 ? 0 : fabs(x[i]));
      }
    }
    euclTrue = 0;
    for (unsigned int i = 0; i < n; ++i) {
      euclTrue += ((ux[i] ? 1. : -1.) - x[i]) * ((ux[i] ? 1. : -1.) - x[i]);
    }
    start = std::chrono::high_resolution_clock::now();

    REQUIRE(decoders.size() == 1);
    calcWeight = decoders[0]->decode(x, y);

    stop = std::chrono::high_resolution_clock::now();
    overall += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    flag = check_decoded_word();
    int syndsum = calculate_syndrome(y);
    euclCalc = 0;
    for (unsigned int i = 0; i < n; ++i) {
      euclCalc += ((y[i] ? 1. : -1.) - x[i]) * ((y[i] ? 1. : -1.) - x[i]);
    }
    INFO("Decoded word: " << arrayToSstream<int>(n, y).str());
    INFO("True weight: " << metric);
    INFO("Calculated weight: " << calcWeight);
    INFO("Syndrome: " << arrayToSstream<int>(n - k, synd).str());
    INFO("True Euclidean metric: " << euclTrue);
    INFO("Resulted Euclidean metric: " << euclCalc);
    // INFO("Check matrix:\n"<<matrixToSstream(n-k, n, g).str());
    if (!((flag || calcWeight < metric || fabs(calcWeight - metric) < EPSILON) &&
        syndsum == 0 &&
        (euclCalc < euclTrue || fabs(euclCalc - euclTrue) < EPSILON))) {
      ++fails;
      CHECK(false);
    }
    if (!((test + 1) % 250)) {
      WARN("Time to decode " << test + 1 << " words: " << overall << "ms");
    }
  }
  WARN("Overall time: " << overall << "ms");
  WARN("Failed tests: " << fails);
}
#endif

void Simulation::log_header() {
  std::cout << "STN,DEV,ITERS,";
  for (auto dec : decoders) {
    std::string name = idToString(dec->get_id());
    std::cout << name << "-FER,";
    std::cout << name << "-CMP,";
    std::cout << name << "-ADD,";
  }
  std::cout << std::endl;
}

void Simulation::log() {
  std::cout << stn << ",";
  std::cout << dev << ",";
  std::cout << test << ",";
  for (unsigned long i = 0; i < decoders.size(); ++i) {
    std::cout << 1.0 * errors_by_decoder[i] / test << "," << 1.0 * cmp_ops_by_decoder[i] / test << ","
              << 1.0 * add_ops_by_decoder[i] / test << ",";
  }
  std::cout << std::endl;
}
std::vector<DecoderID> Simulation::get_decoders() {
  std::vector<DecoderID> res;
  for (auto dec : decoders) {
    res.push_back(dec->get_id());
  }
  return res;
}
unsigned int Simulation::get_max_errors() {
  return max_errors;
}
unsigned int Simulation::get_max_iterations() {
  return max_iterations;
}
unsigned int Simulation::get_code_length() {
  return n;
}
unsigned int Simulation::get_code_dimension() {
  return k;
}
std::string Simulation::get_code_name() {
  return code_name;
}


