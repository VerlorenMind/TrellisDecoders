#include "testutil.h"

TEST_CASE("CONSTRUCTOR: Can read matrix from file")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, n, k);

    REQUIRE(matrix[0] == 0b111);
    REQUIRE(matrix[1] == 0b100);
    REQUIRE(matrix[2] == 0b011);
    REQUIRE(matrix[3] == 0b101);
    REQUIRE(matrix[4] == 0b010);
    REQUIRE(matrix[5] == 0b001);

    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form is achieved")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, 6, 3);

    minspan_form(6, 3, matrix);

    REQUIRE(matrix[0] == 0b001);
    REQUIRE(matrix[1] == 0b011);
    REQUIRE(matrix[2] == 0b010);
    REQUIRE(matrix[3] == 0b111);
    REQUIRE(matrix[4] == 0b110);
    REQUIRE(matrix[5] == 0b100);

    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges can be found")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, n, k);

    unsigned int *ranges = find_ranges(n, k, matrix);

    REQUIRE(ranges[0] == 0);
    REQUIRE(ranges[1] == 5);
    REQUIRE(ranges[2] == 0);
    REQUIRE(ranges[3] == 4);
    REQUIRE(ranges[4] == 0);
    REQUIRE(ranges[5] == 3);

    delete[] matrix;
    delete[] ranges;
}