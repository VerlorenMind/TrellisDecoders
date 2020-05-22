#include <iostream>
#include <chrono>
#include <random>
#include <BSDDecoder.h>
#include <cassert>
#include <Simulation.h>
#include "BeastDecoder.h"

void usage()
{
    std::cout<<"Usage: input_code_file run_file start_stn stop_stn step\ninput_code_file - file containing length, dimension, generator and check matrix of the code\nrun_file - file specifying run configuration\nstart_stn - starting STN value\nstop_stn - final STN value";
}

int main (int argc, char* argv[]){
    if(argc != 6) {
        usage();
    }
    else {
      unsigned int i=1;
      std::ifstream code_input(argv[i++]);
      std::ifstream run_input(argv[i++]);
      double start_stn = strtod(argv[i++], nullptr);
      double stop_stn = strtod(argv[i++], nullptr);
      double step = strtod(argv[i++], nullptr);
      Simulation sim(code_input, run_input);
      sim.log_header();
      for(double stn = start_stn; stn <= stop_stn; stn += step)
      {
        sim.setSTN(stn);
        sim.run();
        sim.log();
      }
    }
    return 0;
}