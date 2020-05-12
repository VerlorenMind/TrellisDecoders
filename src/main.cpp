#include <iostream>
#include <chrono>
#include <random>
#include <BSDDecoder.h>
#include <cassert>
#include "BeastDecoder.h"

void generate_vector(unsigned int size, int* vector)
{
    // Generating random binary vector
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(0,1);
    for(unsigned int i=0; i<size; ++i)
    {
        vector[i] = distribution(generator);
    }
}

void apply_noise(const int* x, double* y, double dev, unsigned int len)
{
    // Applying AWGN with BPSK modulation
    unsigned int seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<double> rv(0.0, dev);
    for(unsigned int i=0; i<len; ++i)
    {
        y[i] = (x[i] ? 1.0 : -1.0) + rv(gen);
    }
}

void usage()
{
    std::cout<<"Usage: beast errors total stn d file file\nerrors - allowed number of errors\ntotal - total decoder runs\nstn - signal-to-noise ratio\nd - delta parameter of the decoder\nfile - file with code parameters, generator and check matrices\n";
}

int main (int argc, char* argv[]){
    if(argc != 6) {
        usage();
    }
    else {
        auto overallStart = std::chrono::high_resolution_clock::now();
        unsigned int maxErrors = (unsigned int) strtol(argv[1], nullptr, 10),
            maxTotal = (unsigned int) strtol(argv[2], nullptr, 10),
            bsdErrorsCount = 0,beastErrorsCount = 0, totalCount = 0;
        double stn = strtod(argv[3], nullptr), delta = strtod(argv[4], nullptr);
        std::ifstream input(argv[5]);
        unsigned int n, k;
        input >> n >> k;
        double dev = sqrt((1 / (((1.0 * k) / n) * pow(10, stn / 10))) / 2);

        int **g = readMatrix(input, n, k);
        int **h = readMatrix(input, n, n - k);

        BeastDecoder beastdec(n, n - k, h);
        BSDDecoder bsddec(n, n - k, h);

        double *x = new double[n];
        int *beasty = new int[n];
        int *bsdy = new int[n];
        int *u = new int[k];
        unsigned long beast_op_add = 0, beast_op_cmp = 0;
        unsigned long bsd_op_add = 0, bsd_op_cmp = 0;
        int *ux = new int[n];
        unsigned int temp;
        bool beastflag, bsdflag;
        std::string word;

        while(totalCount < maxTotal && (bsdErrorsCount < maxErrors && beastErrorsCount < maxErrors)) {
            generate_vector(k, u);
            for (unsigned int i = 0; i < n; ++i) {
                temp = 0;
                for (unsigned int j = 0; j < k; ++j) {
                    temp ^= g[j][i] & u[j];
                }
                ux[i] = (temp ? 1 : 0);
            }
            
            apply_noise(ux, x, dev, n);
            beastdec.decode(x, beasty, delta);
            bsddec.decode(x, bsdy, delta);
            beast_op_add += beastdec.op_add;
            beast_op_cmp += beastdec.op_cmp;
            bsd_op_add += bsddec.op_add;
            bsd_op_cmp += bsddec.op_cmp;
            bsdflag = true;
            beastflag = true;
            for (unsigned int i = 0; i < n; ++i) {
                if (beasty[i] != ux[i]) {
                    beastflag = false;
                }
                if (bsdy[i] != ux[i]) {
                    bsdflag = false;
                }
            }
            if(!bsdflag)
            {
                ++bsdErrorsCount;
            }
            if(!beastflag)
            {
                ++beastErrorsCount;
            }
            ++totalCount;
        }
        delete[] x;
        delete[] beasty;
        delete[] bsdy;
        delete[] u;
        for(unsigned int i=0; i<k; ++i)
        {
            delete[] g[i];
        }
        delete[] g;
        for(unsigned int i=0; i<n-k; ++i)
        {
            delete[] h[i];
        }
        delete[] h;
        delete[] ux;
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout<<stn<<",";
        std::cout<<double(beastErrorsCount)/totalCount<<",";
        std::cout<<double(bsdErrorsCount)/totalCount;
        std::cout<<","<<(((std::chrono::duration<double, std::milli>)(stop - overallStart)).count())/totalCount;
        std::cout<<","<<double(beast_op_add)/totalCount<<","
                      <<double(beast_op_cmp)/totalCount<<","
                      <<double(bsd_op_add)/totalCount<<","
                      <<double(bsd_op_cmp)/totalCount<<","
                      <<std::endl;
    }
}