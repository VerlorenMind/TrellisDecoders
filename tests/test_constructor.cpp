#include <iostream>
#include <stack>
#include <random>
#include <chrono>
#include "catch.hpp"
#include "BeastDecoder.h"
#include "myutil.h"

const double EPSILON = 0.0000001;

TEST_CASE("Can read matrix from file")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int* matrix = readMatrix(filename, 6, 3);

    REQUIRE(matrix[0] == 0b111);
    REQUIRE(matrix[1] == 0b100);
    REQUIRE(matrix[2] == 0b011);
    REQUIRE(matrix[3] == 0b101);
    REQUIRE(matrix[4] == 0b010);
    REQUIRE(matrix[5] == 0b001);

    delete[] matrix;
}

TEST_CASE("Minspan form is achieved")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int* matrix = readMatrix(filename, 6, 3);

    minspan_form(6, 3, matrix);

    REQUIRE(matrix[0] == 0b001);
    REQUIRE(matrix[1] == 0b011);
    REQUIRE(matrix[2] == 0b010);
    REQUIRE(matrix[3] == 0b111);
    REQUIRE(matrix[4] == 0b110);
    REQUIRE(matrix[5] == 0b100);

    delete[] matrix;
}

TEST_CASE("Ranges can be found")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int* matrix = readMatrix(filename, 6, 3);

    unsigned int* ranges = find_ranges(6, 3, matrix);

    REQUIRE(ranges[0] == 0);
    REQUIRE(ranges[1] == 5);
    REQUIRE(ranges[2] == 0);
    REQUIRE(ranges[3] == 4);
    REQUIRE(ranges[4] == 0);
    REQUIRE(ranges[5] == 3);

    delete[] matrix;
    delete[] ranges;
}

TEST_CASE("Can decode zero word")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int* g;
    g = readMatrix(filename, 6, 3);
    BeastDecoder dec(6, 3, filename);
    double* x = new double[6];
    unsigned int* u = new unsigned int[6];
    double delta = 0.5;
    x[0] = -1.5;
    x[1] = -0.5;
    x[2] = -1.7;
    x[3] = -1.9;
    x[4] = -2.0;
    x[5] = -1.0;

    dec.decode(x, u, delta);

    for(unsigned int i=0; i<6; ++i)
    {
        REQUIRE(u[i] == 0);
    }

    delete[] x;
    delete[] u;
    delete[] g;
}

TEST_CASE("Can decode non-zero word")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int* g;
    g = readMatrix(filename, 6, 3);
    BeastDecoder dec(6, 3, filename);
    double* x = new double[6];
    unsigned int* u = new unsigned int[6];
    double delta = 0.5;
    x[0] = 1.5;
    x[1] = 0.5;
    x[2] = -1.7;
    x[3] = -1.9;
    x[4] = 2.0;
    x[5] = 1.0;
    std::string codeword = "000000";

    dec.decode(x, u, delta);

    for(unsigned int j=0; j<6; ++j)
    {
        codeword[j] = u[j] ? '1' : '0';
    }
    INFO("Decoded word: "<<codeword);
    for(unsigned int i=0; i<6; ++i)
    {
        switch(i)
        {
            case 0:
            case 1:
            case 4:
            case 5: {
                INFO("Failed position: "<<i);
                CHECK(u[i] == 1);
                break;
            }
            case 2:
            case 3: {
                INFO("Failed position: "<<i);
                CHECK(u[i] == 0);
                break;
            }
            default:
                break;
        }
    }

    delete[] x;
    delete[] u;
}

void generate_vector(unsigned int size, unsigned int* vector)
{
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(0,1);
    for(unsigned int i=0; i<size; ++i)
    {
        vector[i] = distribution(generator);
    }
}

TEST_CASE("Can decode series of random words")
{
    double stn = 100, delta = 0.5;
    unsigned int n = 6, k=3;
    double dev = sqrt((1/(((1.0 * k) / n)*pow(10, stn / 10)))/2);
    unsigned int tests = 1000;

    unsigned int seed = (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<double> rv(0.0, (double) sqrt(dev));

    std::ifstream input("../tests/test_matrix");
    unsigned int* g = readMatrix(input, n, k);
    BeastDecoder dec(n, k, input);
    double* x = new double[n];
    unsigned int* y = new unsigned int[n];
    unsigned int* u = new unsigned int[k];
    double trueWeight, calcWeight, noise;
    unsigned int temp;
    bool flag;
    std::string word;


    for(unsigned int test = 0; test < tests; ++test)
    {
        trueWeight = 0;
        generate_vector(k, u);
        INFO("Informational word: " << array_to_sstream<unsigned int>(k, u).str());
        for(unsigned int i=0; i<n; ++i)
        {
            temp = 0;
            for(unsigned int j=0; j<k; ++j)
            {
                temp ^= ((g[i] & (u[j] << j))>>j);
            }
            noise = rv(gen);
            x[i] = (temp ? 1. : -1.) + noise;
            trueWeight += fabs(noise);
        }
        INFO("Coded word with noise: " << array_to_sstream<double>(n, x).str());
        calcWeight = dec.decode(x, y, delta);
        flag = true;
        for(unsigned int i=0; i<n; ++i)
        {
            if(y[i] != u[i])
            {
                flag = false;
                break;
            }
        }
        INFO("Decoded word: " << array_to_sstream<unsigned int>(n, y).str());
        INFO("True weight: " << trueWeight);
        INFO("Calculated weight: " << calcWeight);
        CHECK((flag || fabs(calcWeight - trueWeight) < EPSILON));
    }

    delete[] x;
    delete[] y;
    delete[] u;
    delete[] g;
}