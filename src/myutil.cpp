#include <iostream>
#include "myutil.h"


int** readMatrix(std::ifstream& input, unsigned int n, unsigned int k)
{
    // Reading matrix from a file and storing it's columns as uints
    int** matrix = new int*[k];
    for(unsigned int i=0; i<k; ++i)
    {
        matrix[i] = new int[n];
    }
    uint64_t bit;
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            input >> matrix[i][j];
        }
    }
    return matrix;
}

void minspan_form(unsigned int n, unsigned int k, int** a) {
    unsigned int num; // number of rows with 1s in i'th column
    unsigned int fixed_rows = 0; // number of rows on top of the matrix that were transformed
    unsigned int* rows = new unsigned int[k]; // rows with 1s in i'th column
    uint64_t temp;
    // Left side
    for(unsigned int i=0; i<k; ++i)
    {
        num = 0;
        // Finding all rows with 1s in i'th column
        for(unsigned int j=fixed_rows; j<k; ++j)
        {
            if(a[j][i])
            {
                rows[num++] = j;
            }
        }
        if(num == 0)
        {
            continue; // should never happen with non-zero determinant TODO: handle this exception
        }
        else
        {
            // If the first row after fixed rows does not contain 1, swap it with the first row that does TODO: rewrite this part
            if(rows[0] != fixed_rows)
            {
                for(unsigned int l=0; l<n; ++l)
                {
                   // temp = a[l] & (uint64_t(1) << rows[0]);
                   // a[l] ^= (((a[l] & (uint64_t(1) << i)) << (rows[0]-i)) ^ (a[l] & (uint64_t(1) << rows[0])));
                   // a[l] ^= (temp >> (rows[0]-i)) ^ (a[l] & (uint64_t(1)<<i));
                   std::swap(a[rows[0]][l], a[i][l]);
                }
                rows[0] = i;
            }

            ++fixed_rows;

            // Adding fixed row to the unfixed rest of the matrix
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    // a[j] ^= (a[j] & (uint64_t(1)<<rows[0]))<<(rows[l]-rows[0]);
                    a[rows[l]][j] ^= a[rows[0]][j];
                }
            }
        }
    }
    // Right side
    // Same stuff as above, but with different indices and without swapping rows
    fixed_rows = 0;
    int *fixed_nums = new int[n-k];
    std::string tempmatr;
    int i = n-1;
    while(fixed_rows < k)
    {
        tempmatr = matrix_to_sstream(k, n, a).str();
        num = 0;
        for(int j=k-1; j>=0; --j)
        {
            bool flag = true;
            for(unsigned int l=0; l<fixed_rows; ++l)
            {
                if(fixed_nums[l] == j)
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                if (a[j][i])
                {
                    rows[num++] = (unsigned int) j;
                }
            }
        }
        --i;
        if(num == 0)
        {
            continue;
        }
        else
        {
            fixed_nums[fixed_rows++] = rows[0];
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[rows[l]][j] ^= a[rows[0]][j];
                }
            }
        }
    }
    delete [] rows;
}

unsigned int* find_ranges(unsigned int n, unsigned int k, int** a)
{
    // Finding active ranges for each row (range between the first non-zero element and the last one)
    unsigned int* ranges = new unsigned int[2*k];
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            if(a[i][j])
            {
                ranges[2*i] = j;
                break;
            }
        }
        for(int j=n-1; j>=0; --j) {
            if (a[i][j]) {
                ranges[2*i+1] = (unsigned int) j;
                break;
            }
        }
    }
    return ranges;
}


