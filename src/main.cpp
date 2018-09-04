#include <iostream>
#include <chrono>
#include <random>
#include <myutil.h>
#include "BeastDecoder.h"

void generate_vector(unsigned int size, unsigned int* vector)
{
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(0,1);
    for(unsigned int i=0; i<size; ++i)
    {
        vector[i] = (unsigned int) distribution(generator);
    }
}

void apply_noise(const int* x, double* y, double dev, unsigned int len)
{
    unsigned int seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<double> rv(0.0, sqrt(dev));
    for(unsigned int i=0; i<len; ++i)
    {
        y[i] = (x[i] ? 1.0 : -1.0) + rv(gen);
    }
}

void usage()
{
    std::cout<<"Usage: beast errors total stn d file\nerrors - allowed number of errors\ntotal - total decoder runs\nstn - signal-to-noise ratio\nd - delta parameter of the decoder\nfile - file with code parameters, generator and check matrices\n";
}

int main (int argc, char* argv[]){
    if(argc != 6) {
        usage();
    }
    else {
        auto overallStart = std::chrono::high_resolution_clock::now();
        unsigned int maxErrors = (unsigned int) strtol(argv[1], nullptr, 10),
            maxTotal = (unsigned int) strtol(argv[2], nullptr, 10),
            errorsCount, totalCount;
        double stn = strtod(argv[3], nullptr), delta = strtod(argv[4], nullptr);
        std::ifstream input(argv[5]);
        unsigned int n, k;
        input >> n >> k;
        double dev = sqrt((1 / (((1.0 * k) / n) * pow(10, stn / 10))) / 2);

        uint64_t *g = readMatrix(input, n, k);

        BeastDecoder dec(n, n - k, input);

        double *x = new double[n];
        unsigned int *y = new unsigned int[n];
        unsigned int *u = new unsigned int[k];
        int *ux = new int[n];
        double trueWeight, calcWeight, noise;
        unsigned int temp;
        bool flag;
        std::string word;

        while(totalCount < maxTotal && errorsCount < maxErrors) {
            trueWeight = 0;
            generate_vector(k, u);
            for (unsigned int i = 0; i < n; ++i) {
                temp = 0;
                for (unsigned int j = 0; j < k; ++j) {
                    temp ^= ((g[i] & (u[j] << j)) >> j);
                }
                ux[i] = (temp ? 1 : 0);
            }
            apply_noise(ux, x, dev, n);
            calcWeight = dec.decode(x, y, delta);
            flag = true;
            for (unsigned int i = 0; i < n; ++i) {
                if (y[i] != ux[i]) {
                    flag = false;
                    break;
                }
            }
            if(!flag)
            {
                ++errorsCount;
            }
            ++totalCount;
        }
        delete[] x;
        delete[] y;
        delete[] u;
        delete[] g;
        delete[] ux;
        std::cout<<errorsCount/totalCount<<std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        std::cout<<"Execution time: "<<((std::chrono::duration<double, std::milli>)(stop - overallStart)).count()<<"ms\n";
    }
}