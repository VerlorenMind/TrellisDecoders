#ifndef BEAST_SIMULATION_H
#define BEAST_SIMULATION_H

#include <vector>
#include <random>
#include "CatchWrap.h"
#include "DecoderID.h"
#include "SoftDecoder.h"

class Simulation {
 private:
  unsigned int n, k, test;
  unsigned int max_iterations, max_errors;
  int **g, **h;
  double dev, stn;
  double* x;
  int* y;
  int* u;
  int* ux;
  int *synd;
  std::string code_name;
  std::vector<SoftDecoder*> decoders;
  std::default_random_engine gen;
  void generate_information_vector();
  void encode();
  void apply_noise();
  int calculate_syndrome(const int* vec);
  bool check_decoded_word();
 public:
  std::vector<long> errors_by_decoder;
  std::vector<long> cmp_ops_by_decoder;
  std::vector<long> add_ops_by_decoder;
  explicit Simulation(std::ifstream &input_file, unsigned int max_errors = 0,
                      unsigned int max_iterations = 0, unsigned int seed = 123);
  Simulation(std::ifstream &input_file, std::ifstream &run_config);
  Simulation(std::ifstream &input_file, double stn, SoftDecoder* decoder, unsigned int max_errors,
             unsigned int max_iterations, unsigned int seed = 123);
  ~Simulation();
  void run();
  SoftDecoder* get_decoder(DecoderID id);
  void add_decoder(DecoderID id);
  void setSTN(double stn);
#ifdef CATCH_TESTING
  void test_run();
#endif
  void log();
  void log_header();
  std::vector<DecoderID> get_decoders();
  unsigned int get_max_errors();
  unsigned int get_max_iterations();
  unsigned int get_code_length();
  unsigned int get_code_dimension();
  std::string get_code_name();
};

void generate_vector(unsigned int size, int* vector);
#endif //BEAST_SIMULATION_H
