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


/*TEST_CASE("Simple trellis can be constructed")
{
    bool genmat_c[3][6] = {
            {1, 1, 1, 0, 0, 0},
            {0, 1, 1, 1, 1, 0},
            {0, 0, 1, 0, 1, 1}
    };
    bool** genmat = new bool*[3];
    genmat[0] = (bool*)&(genmat_c[0]);
    genmat[1] = (bool*)&(genmat_c[1]);
    genmat[2] = (bool*)&(genmat_c[2]);
    Node** result;
    std::deque<Node> outputDeque;

    result = construct_trellis(6, 3, genmat);

    REQUIRE(result[0][0].number == 0);
    REQUIRE(result[6][0].number == 0);
    printTree(result[0][0]);

    delete [] genmat;
    for(unsigned int i=0; i<7; ++i)
    {
        delete [] result[i];
    }
    delete [] result;

}*/