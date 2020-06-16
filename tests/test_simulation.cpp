#include <fstream>
#include <BeastDecoder.h>
#include "catch.hpp"
#include "Simulation.h"

const double EPSILON = 0.0000001;

TEST_CASE("SIM: Can create empty simulation") {
  std::ifstream in("../tests/test_matrix");
  Simulation sim(in, 10, 10, 0);

  REQUIRE(sim.get_code_name() == "Test matrix");
  REQUIRE(sim.get_code_length() == 6);
  REQUIRE(sim.get_code_dimension() == 3);
  REQUIRE(sim.get_max_errors() == 10);
  REQUIRE(sim.get_max_iterations() == 10);
  REQUIRE(sim.get_decoders().empty());
}

TEST_CASE("SIM: Can add decoder to the simulation") {
  std::ifstream in("../tests/test_matrix");
  Simulation sim(in, 10, 10, 0);
  sim.add_bsd_decoder();
  SoftDecoder *dec = sim.get_decoder(DecoderID::BSD);

  REQUIRE(dec != nullptr);
  REQUIRE(sim.get_decoders().size() == 1);
}

TEST_CASE("SIM: Can read spec file") {
  std::ifstream code_in("../tests/test_matrix");
  std::ifstream run_in("../data/run_all");
  Simulation sim(code_in, run_in);
  SoftDecoder *beast = sim.get_decoder(DecoderID::BEAST);
  SoftDecoder *bsd = sim.get_decoder(DecoderID::BSD);

  REQUIRE(sim.get_max_errors() == 1000);
  REQUIRE(sim.get_max_iterations() == 1000000);
  REQUIRE(sim.get_decoders().size() == 2);
  REQUIRE(beast != nullptr);
  REQUIRE(bsd != nullptr);
  double temp = fabs((dynamic_cast<BeastDecoder *>(beast)->get_delta() - 0.5));
  REQUIRE(temp < EPSILON);
}

TEST_CASE("SIM: Can launch a simple run") {
  std::ifstream in("../tests/test_matrix");
  Simulation sim(in, 10, 10, 0);
  sim.add_bsd_decoder();
  sim.setSTN(100);

  sim.run();

  REQUIRE(sim.errors_by_decoder[0] == 0);
  REQUIRE(sim.cmp_ops_by_decoder[0] > 0);
  REQUIRE(sim.add_ops_by_decoder[0] > 0);
}