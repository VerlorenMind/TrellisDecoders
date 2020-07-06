#ifndef BEAST_SIMULATION_H
#define BEAST_SIMULATION_H

#include <vector>
#include <random>
#include "CatchWrap.h"
#include "DecoderID.h"
#include "SoftDecoder.h"

#ifdef CATCH_TESTING
struct TestCase {
  int *ux, *y;
  double *x;
};
#endif
class Simulation {
 private:
  unsigned int n, k, test;
  unsigned int max_iterations, max_errors;
  int **g, **h;
  double dev, stn;
  double *x;
  int *y;
  int *u;
  int *ux;
  int *synd;
  std::string code_name;
  std::vector<SoftDecoder *> decoders;
  std::default_random_engine gen;
  void generate_information_vector();
  void encode();
  void apply_noise();
  int calculate_syndrome(const int *vec);
  bool check_decoded_word(int dec_num);
  bool check_decoded_word();
 public:
#ifdef CATCH_TESTING
 std::vector<TestCase> failed_cases;
#endif
  std::vector<long> errors_by_decoder;
  std::vector<long> iters_by_decoder;
  std::vector<long> cmp_ops_by_decoder;
  std::vector<long> add_ops_by_decoder;
  std::vector<long> bit_errors_by_decoder;
  std::vector<long> algo_iters_by_decoder;
  explicit Simulation(std::ifstream &input_file, unsigned int max_errors = 0,
                      unsigned int max_iterations = 0, unsigned int seed = 123);
  Simulation(std::ifstream &input_file, std::ifstream &run_config);
  Simulation(std::ifstream &input_file, double stn, SoftDecoder *decoder, unsigned int max_errors,
             unsigned int max_iterations, unsigned int seed = 123);
  ~Simulation();
  void run();
  SoftDecoder *get_decoder(DecoderID id);
  // void add_decoder(DecoderID id);
  void add_beast_decoder(double delta);
  void add_bsd_decoder();
  void add_viterbi_decoder();
  void add_ktkl_decoder(int w, int buf_size);
  void add_ordered_statistics_decoder(int order);
  void setSTN(double stn);
#ifdef CATCH_TESTING
  void test_run(bool silent=false);
  void run_failed_cases(bool silent=false);
#endif
  void log();
  void log_header();
  void clear_decoders();
  std::vector<DecoderID> get_decoders();
  unsigned int get_max_errors();
  unsigned int get_max_iterations();
  unsigned int get_code_length();
  unsigned int get_code_dimension();
  std::string get_code_name();
};

void generate_vector(unsigned int size, int *vector);
#endif //BEAST_SIMULATION_H
