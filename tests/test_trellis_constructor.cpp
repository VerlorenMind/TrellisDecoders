#include "testutil.h"

TEST_CASE("CONSTRUCTOR: Can read test matrix from file")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);
    int answer[3][6] = {
            {1, 0, 1, 1, 0, 1},
            {1, 0, 1, 0, 1, 0},
            {1, 1, 0, 1, 0, 0},
    };
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            REQUIRE(answer[i][j] == matrix[i][j]);
        }
    }
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form for test matrix is achieved")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);
    int answer[3][6] = {
            {1, 1, 0, 1, 0, 0},
            {0, 1, 1, 1, 1, 0},
            {0, 0, 0, 1, 1, 1},
    };
    minspan_form(n, k, matrix);

    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            REQUIRE(answer[i][j] == matrix[i][j]);
        }
    }
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for  test matrix can be found")
{
    std::ifstream filename("../tests/test_matrix");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);

    unsigned int *ranges = find_ranges(n, k, matrix);

    REQUIRE(ranges[0] == 0);
    REQUIRE(ranges[1] == 5);
    REQUIRE(ranges[2] == 0);
    REQUIRE(ranges[3] == 4);
    REQUIRE(ranges[4] == 0);
    REQUIRE(ranges[5] == 3);

    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] ranges;
}

TEST_CASE("CONSTRUCTOR: Can read RM(16,5) matrix from file")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);

    int answer[5][16] = {
            {1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1},
            {0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1},
            {0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1},
            {0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1},
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
    };
    for (unsigned int i = 0; i < k; ++i)
    {
        for (unsigned int j = 0; j < n; ++j)
        {
            REQUIRE(answer[i][j] == matrix[i][j]);
        }
    }
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Minspan form for RM(16,5) matrix is achieved")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);

    int answer[5][16] = {
            {1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0},
            {0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0},
            {0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0},
            {0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
    };
    minspan_form(n, k, matrix);

    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            REQUIRE(answer[i][j] == matrix[i][j]);
        }
    }
    // WARN("Resulted matrix:\n" << matrix_to_sstream(k, n, matrix).str());
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

TEST_CASE("CONSTRUCTOR: Ranges for RM(16, 5) matrix can be found")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);
    minspan_form(n, k, matrix);

    unsigned int *ranges = find_ranges(n, k, matrix);

    // WARN("Resulted matrix:\n" << matrix_to_sstream(k, n, matrix).str());
    REQUIRE(ranges[0] == 0);
    REQUIRE(ranges[1] == 8);
    REQUIRE(ranges[2] == 1);
    REQUIRE(ranges[3] == 14);
    REQUIRE(ranges[4] == 2);
    REQUIRE(ranges[5] == 13);
    REQUIRE(ranges[6] == 3);
    REQUIRE(ranges[7] == 11);
    REQUIRE(ranges[8] == 4);
    REQUIRE(ranges[9] == 15);

    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] ranges;
}

TEST_CASE("CONSTRUCTOR: Minspan form for RM(16,5) check matrix is achieved")
{
    std::ifstream filename("../data/reed-muller-16-1");
    unsigned int n, k;
    filename >> n >> k;
    int **matrix = readMatrix(filename, n, k);
    // WARN("Read matrix:\n" << matrix_to_sstream<int>(k, n, matrix).str());
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
    matrix = readMatrix(filename, n, n-k);
    // WARN("Read matrix:\n" << matrix_to_sstream<int>(n-k, n, matrix).str());

    int answer[11][16] = {
            {1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1},
    };
    minspan_form(n, n-k, matrix);

    // WARN("Resulted matrix:\n" << matrix_to_sstream(n-k, n, matrix).str());
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            REQUIRE(answer[i][j] == matrix[i][j]);
        }
    }
    for (unsigned int i = 0; i < k; ++i)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}
