#define TESTING

#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include <cmath>
#include "myutil.h"
#include "BeastDecoder.h"

unsigned int* readMatrix(const char* filename, unsigned int n, unsigned int k)
{
    // Reading matrix from a file and storing it's columns as uints
    std::ifstream input(filename);
    unsigned int* matrix = new unsigned int[n];
    for(unsigned int i=0; i<n; ++i)
    {
        matrix[i] = 0;
    }
    unsigned int bit;
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            input >> bit;
            matrix[j] ^= (bit << i);
        }
    }
    return matrix;
}

void minspan_form(unsigned int n, unsigned int k, unsigned int* a) {
    unsigned int num; // number of rows with 1s in i'th column
    unsigned int fixed_rows = 0; // number of rows on top of the matrix that were transformed
    unsigned int* rows = new unsigned int[k]; // rows with 1s in i'th column
    // Left side
    for(unsigned int i=0; i<k; ++i)
    {
        num = 0;
        // Finding all rows with 1s in i'th column
        for(unsigned int j=fixed_rows; j<k; ++j)
        {
            if(a[i] & (1<<j))
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
            // If the first row after fixed rows does not contain 1, swap it with the first row that does
            if(rows[0] != fixed_rows)
            {
                unsigned int temp;
                for(unsigned int l=0; l<n; ++l)
                {
                    // Swapping bits
                    temp = a[l] & (1 << rows[0]);
                    a[l] ^= (((a[l] & (1 << i)) << (rows[0]-i)) ^ (a[l] & (1 << rows[0])));
                    a[l] ^= (temp >> (rows[0]-i)) ^ (a[l] & (1<<i));
                }
                rows[0] = i;
            }

            ++fixed_rows;

            // Adding fixed row to the unfixed rest of the matrix
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[j] ^= (a[j] & (1<<rows[0]))<<(rows[l]-rows[0]);
                }
            }
        }
    }
    // Right side
    // Same stuff as above, but with different indices and without swapping rows
    fixed_rows = 0;
    for(unsigned int i=n-1; i>=n-k; --i)
    {
        num = 0;
        for(int j=k-fixed_rows-1; j>=0; --j)
        {
            if(a[i] & (1<<j))
            {
                rows[num++] = j;
            }
        }
        if(num == 0)
        {
            continue;
        }
        else
        {
            ++fixed_rows;
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[j] ^= (a[j] & (1<<rows[0]))>>(rows[0]-rows[l]);
                }
            }
        }
    }
    delete [] rows;
}

unsigned int* find_ranges(unsigned int n, unsigned int k, unsigned int* a)
{
    // Finding active ranges for each row (range between the first non-zero element and the last one)
    unsigned int* ranges = new unsigned int[2*k];
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            if(a[j] & (1<<i))
            {
                ranges[2*i] = j;
                break;
            }
        }
        for(unsigned int j=n-1; j>=0; --j) {
            if (a[j] & (1<<i)) {
                ranges[2*i+1] = j;
                break;
            }
        }
    }
    return ranges;
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, const char* matrix_file)
{
    // Transforming matrix to minspan form and finding ative ranges of it's rows
    this->n = n;
    this->k = k;
    h = readMatrix(matrix_file, n, k);
    minspan_form(n, k, h);
    ranges = find_ranges(n, k, h);
}

BeastDecoder::~BeastDecoder()
{
    delete[] h;
    delete[] ranges;
}

inline double metric(double x, int y)
{
    return fabs(y - x);
}

void BeastDecoder::decode(double *x, unsigned int *u, double delta)
{
    std::list<Node> fwdTree, bkwTree; // forward and backward trees
    // Initializing starting nodes
    Node temp;
    temp.number = 0;
    temp.layer = 0;
    temp.metric = 0;
    temp.path.resize(0);
    temp.path0 = true;
    temp.path1 = true;
    fwdTree.push_back(temp);
    temp.number = 0;
    temp.layer = n;
    temp.metric = 0;
    temp.path.resize(0);
    temp.path0 = true;
    temp.path1 = true;
    bkwTree.push_back(temp);
    unsigned int layerMask; // a mask that denotes zero bit-positions on a layer
    double metricBound = 0; // target metric bound
    double min_metric = -1;
    std::vector<unsigned int> min_candidate;
    while(min_metric == -1) {
        metricBound += delta;
        // Growing forward tree
        for (auto iter = fwdTree.begin(); iter != fwdTree.end(); ++iter) {
            // Getting layerMask for current layer
            layerMask = 0;
            for (unsigned int i = 0; i < k; ++i) {
                if (iter->layer < ranges[2 * i] || iter->layer >= ranges[2 * i + 1]) {
                    layerMask ^= (1 << i);
                }
            }
            // If transition with input 0 is possible, add node to the list
            if (iter->path0 && !(iter->number & layerMask) && (iter->metric < metricBound)) {
                temp.number = iter->number;
                temp.layer = iter->layer + 1;
                temp.metric = iter->metric + metric(x[iter->layer], -1);
                temp.path = iter->path;
                temp.path.push_back(0);
                // Inserting new node in sorted list
                for (auto iterSort = iter;; ++iterSort) {
                    if (iterSort == fwdTree.end() ||
                        (iterSort->layer == temp.layer && iterSort->number > temp.number) ||
                        (iterSort->layer > temp.layer)) {
                        fwdTree.insert(iterSort, temp);
                        break;
                    }
                }
                iter->path0 = false;
            }
            // Same for input 1
            if (iter->path1 && !((iter->number ^ h[iter->layer]) & layerMask) && iter->metric < metricBound) {
                temp.number = iter->number ^ h[iter->layer];
                temp.layer = iter->layer + 1;
                temp.metric = iter->metric + metric(x[iter->layer], 1);
                temp.path = iter->path;
                temp.path.push_back(1);
                for (auto iterSort = iter;; ++iterSort) {
                    if (iterSort == fwdTree.end() ||
                        (iterSort->layer == temp.layer && iterSort->number > temp.number) ||
                        (iterSort->layer > temp.layer)) {
                        fwdTree.insert(iterSort, temp);
                        break;
                    }
                }
                iter->path1 = false;
            }
        }
        // Growing backward tree
        for (auto iter = bkwTree.rbegin(); iter != bkwTree.rend(); ++iter) {
            layerMask = 0;
            for (unsigned int i = 0; i < k; ++i) {
                if ((iter->layer-1) <= ranges[2 * i] || (iter->layer-1) > ranges[2 * i + 1]) {
                    layerMask ^= (1 << i);
                }
            }
            // The difference between forward tree is in third condition
            if (iter->path0 && !(iter->number & layerMask) && (iter->metric + metric(x[iter->layer-1], -1)) < metricBound) {
                temp.number = iter->number;
                temp.layer = iter->layer - 1;
                temp.metric = iter->metric + metric(x[iter->layer-1], -1);
                temp.path = iter->path;
                temp.path.insert(temp.path.begin(), 0);
                for (auto iterSort = bkwTree.begin();; ++iterSort) {
                    if (iterSort == bkwTree.end() ||
                        (iterSort->layer == temp.layer && iterSort->number > temp.number) ||
                        (iterSort->layer > temp.layer)) {
                        bkwTree.insert(iterSort, temp);
                        break;
                    }
                }
                iter->path0 = false;
            }
            if (iter->path1 && !((iter->number ^ h[iter->layer-1]) & layerMask) &&
                (iter->metric + metric(x[iter->layer-1], 1)) < metricBound) {
                temp.number = iter->number ^ h[iter->layer];
                temp.layer = iter->layer - 1;
                temp.metric = iter->metric + metric(x[iter->layer-1], 1);
                temp.path = iter->path;
                temp.path.insert(temp.path.begin(), 1);
                for (auto iterSort = bkwTree.begin();; ++iterSort) {
                    if (iterSort == bkwTree.end() ||
                        (iterSort->layer == temp.layer && iterSort->number > temp.number) ||
                        (iterSort->layer > temp.layer)) {
                        bkwTree.insert(iterSort, temp);
                        break;
                    }
                }
                iter->path1 = false;
            }
        }
        // Looking for matches in sorted lists
        auto fwdIter = fwdTree.begin();
        auto bkwIter = bkwTree.begin();
        double tempMetric;
        while(fwdIter != fwdTree.end() && bkwIter != bkwTree.end())
        {
            if(fwdIter->layer == bkwIter->layer && fwdIter->number == bkwIter->number)
            {
                tempMetric = fwdIter->metric + bkwIter->metric;
                if(min_metric == -1 || tempMetric < min_metric) {
                    // Found a match, storing it
                    min_metric = tempMetric;
                    min_candidate.resize(fwdIter->path.size() + bkwIter->path.size());
                    min_candidate.insert(min_candidate.end(), fwdIter->path.begin(), fwdIter->path.end());
                    min_candidate.insert(min_candidate.end(), bkwIter->path.begin(), bkwIter->path.end());
                }
                ++fwdIter;
                ++bkwIter;
            }
            else if((fwdIter->layer < bkwIter->layer) ||
                    (fwdIter->layer == bkwIter->layer && fwdIter->number < bkwIter->number))
            {
                ++fwdIter;
            }
            else
            {
                ++bkwIter;
            }
        }
    }
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = min_candidate[i];
    }
}
