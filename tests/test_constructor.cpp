#include <iostream>
#include <stack>
#include "catch.hpp"
#include "BeastDecoder.h"
#include "myutil.h"


TEST_CASE("Can read matrix from file")
{
    unsigned int* matrix = readMatrix("../tests/test_matrix", 6, 3);

    REQUIRE(matrix[0] == 1);
    REQUIRE(matrix[1] == 5);
    REQUIRE(matrix[2] == 2);
    REQUIRE(matrix[3] == 4);
    REQUIRE(matrix[4] == 3);
    REQUIRE(matrix[5] == 7);

    delete[] matrix;
}

TEST_CASE("Minspan form is achieved")
{
    unsigned int* matrix = readMatrix("../tests/test_matrix", 6, 3);

    minspan_form(6, 3, matrix);

    REQUIRE(matrix[0] == 1);
    REQUIRE(matrix[1] == 3);
    REQUIRE(matrix[2] == 7);
    REQUIRE(matrix[3] == 2);
    REQUIRE(matrix[4] == 6);
    REQUIRE(matrix[5] == 4);

    delete[] matrix;
}

TEST_CASE("Ranges can be found")
{
    unsigned int* matrix = readMatrix("../tests/test_matrix", 6, 3);

    unsigned int* ranges = find_ranges(6, 3, matrix);

    REQUIRE(ranges[0] == 0);
    REQUIRE(ranges[1] == 5);
    REQUIRE(ranges[2] == 2);
    REQUIRE(ranges[3] == 5);
    REQUIRE(ranges[4] == 1);
    REQUIRE(ranges[5] == 5);

    delete[] matrix;
    delete[] ranges;
}

TEST_CASE("Can decode simple word")
{
    BeastDecoder dec(6, 3, "../tests/test_matrix");
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
        REQUIRE(u[0] == 0);
    }

    delete[] x;
    delete[] u;
}