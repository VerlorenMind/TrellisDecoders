#include "testutil.h"

const double EPSILON = 0.0000001;

TEST_CASE("BEAST: Can decode zero word")
{
    std::ifstream filename("../tests/test_matrix");
    uint64_t* g;
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

TEST_CASE("BEAST: Can decode non-zero word")
{
    // notice - with this data it should fail
    std::ifstream filename("../tests/test_matrix");
    uint64_t* g;
    unsigned int n, k;
    filename >> n >> k;
    g = readMatrix(filename, n, k);
    BeastDecoder dec(n, k, filename);
    double* x = new double[6];
    unsigned int* y = new unsigned int[6];
    unsigned int* u = new unsigned int[6];
    double delta = 0.5;
    u[0] = 0;
    u[1] = 1;
    u[2] = 1;
    u[3] = 0;
    u[4] = 0;
    u[5] = 1;
    x[0] = 0.503536;
    x[1] = 0.423184;
    x[2] = 1.87715;
    x[3] = -0.107348;
    x[4] = -1.7087;
    x[5] = -1.10187;
    // True weight: 1.6324

    std::string codeword = "000000";

    dec.decode(x, y, delta);

    for(unsigned int j=0; j<6; ++j)
    {
        codeword[j] = y[j] ? '1' : '0';
    }
    INFO("Decoded word: "<<codeword);
    for(unsigned int i=0; i<6; ++i)
    {
        CHECK(y[i] == u[i]);
    }

    delete[] g;
    delete[] x;
    delete[] u;
    delete[] y;
}

TEST_CASE("BEAST: Can decode series of random words with minimal noise")
{
    test_decoder(1000, 100, 0.5, "../tests/test_matrix", "beast");
}

TEST_CASE("BEAST: Can decode series of random words with some noise")
{
    test_decoder(1000, 1, 0.5, "../tests/test_matrix", "beast");
}

TEST_CASE("BEAST: Can decode BCH(31, 16, 7)")
{
    test_decoder(1000, 1, 0.5, "../data/bch-31-16-7", "beast");
}

TEST_CASE("BEAST: Can decode BCH(31, 21, 5)")
{
    test_decoder(1000, 1, 0.5, "../data/bch-31-21-5", "beast");
}

TEST_CASE("BEAST: Can decode BCH(63, 7, 57)")
{
    test_decoder(1000, 5, 0.5, "../data/bch-63-7-57", "beast");
}

TEST_CASE("BEAST: Can decode BCH(63, 16, 23)")
{
    test_decoder(1000, 1, 0.5, "../data/bch-63-16-23", "beast");
}

TEST_CASE("BEAST: Can decode BCH(63, 30, 13)")
{
    test_decoder(1000, 1, 0.5, "../data/bch-63-30-13", "beast");
}

TEST_CASE("BEAST: Can decode BCH(63, 39, 9)")
{
    test_decoder(1000, 1, 0.5, "../data/bch-63-39-9", "beast");
}