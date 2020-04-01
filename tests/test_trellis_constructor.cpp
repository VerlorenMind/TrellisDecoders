#include "testutil.h"

TEST_CASE("CONSTRUCTOR: Can read test matrix from file")
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

TEST_CASE("CONSTRUCTOR: Minspan form for test matrix is achieved")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, n, k);

    minspan_form(n, k, matrix);

    REQUIRE(matrix[0] == 0b001);
    REQUIRE(matrix[1] == 0b011);
    REQUIRE(matrix[2] == 0b010);
    REQUIRE(matrix[3] == 0b111);
    REQUIRE(matrix[4] == 0b110);
    REQUIRE(matrix[5] == 0b100);

    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for  test matrix can be found")
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

TEST_CASE("CONSTRUCTOR: Can read RM(16,5) matrix from file")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, n, k);

    uint64_t answer[16] =
            {1, 10, 100, 1000,
             10000, 1011, 1101, 1110,
             111, 10011, 10101, 10110,
             11001, 11010, 11100, 11111};
    for(unsigned int i=0; i<n; ++i)
    {
        REQUIRE(answer[i] == matrix[i]);
    }
    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form for RM(16,5) matrix is achieved")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    uint64_t *matrix = readMatrix(filename, n, k);

    minspan_form(n, k, matrix);

    REQUIRE(matrix[0] == 0b001);
    REQUIRE(matrix[1] == 0b011);
    REQUIRE(matrix[2] == 0b010);
    REQUIRE(matrix[3] == 0b111);
    REQUIRE(matrix[4] == 0b110);
    REQUIRE(matrix[5] == 0b100);

    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for RM(16, 5) matrix can be found")
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
