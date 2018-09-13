#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include <cmath>
#include <set>
#include "myutil.h"
#include "BeastDecoder.h"

uint64_t* readMatrix(std::ifstream& input, unsigned int n, unsigned int k)
{
    // Reading matrix from a file and storing it's columns as uints
    uint64_t* matrix = new uint64_t[n];
    for(unsigned int i=0; i<n; ++i)
    {
        matrix[i] = 0;
    }
    uint64_t bit;
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

void minspan_form(unsigned int n, unsigned int k, uint64_t* a) {
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
            // If the first row after fixed rows does not contain 1, sum it with the first row that does
            if(rows[0] != fixed_rows)
            {
                for(unsigned int l=0; l<n; ++l)
                {
                    // a[l] ^= (a[l] & (1 << rows[0])) >> (rows[0] - fixed_rows);
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

unsigned int* find_ranges(unsigned int n, unsigned int k, uint64_t* a)
{
    // Finding active ranges for each row (range between the first non-zero element and the last one)
    unsigned int* ranges = new unsigned int[2*k];
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            if(a[j] & (uint64_t(1)<<i))
            {
                ranges[2*i] = j;
                break;
            }
        }
        for(int j=n-1; j>=0; --j) {
            if (a[j] & (uint64_t(1)<<i)) {
                ranges[2*i+1] = j;
                break;
            }
        }
    }
    return ranges;
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, std::ifstream& filename)
{
    // Transforming matrix to minspan form and finding active ranges of it's rows
    this->n = n;
    this->k = k;
    h = readMatrix(filename, n, k);
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

void insertNode(Node& node, std::set<Node, NodeCompare>& tree)
{
    // Check if this node is already in the tree
    auto tempfind = tree.find(node);
    // If it isn't, add it
    if (tempfind == tree.end()) {
        tree.insert(node);
    } else {
        // If it is, try to minimize it's metric
        if (tempfind->metric > node.metric) {
            tempfind->metric = node.metric;
            tempfind->path = node.path;
        }
    }
}

double BeastDecoder::decode(double *x, unsigned int *u, double delta)
{
    std::set<Node, NodeCompare> *fwdTree, *bkwTree; // forward and backward trees
    std::list<Node>* nodeBuffer = new std::list<Node>[n+1];
    // Initializing starting nodes
    fwdTree = new std::set<Node, NodeCompare>[n+1];
    bkwTree = new std::set<Node, NodeCompare>[n+1];
    Node temp{};
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.path0 = true;
    temp.path1 = true;
    fwdTree[0].insert(temp);
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.path0 = true;
    temp.path1 = true;
    bkwTree[n].insert(temp);
    unsigned int layerMask; // a mask that denotes zero bit-positions on a layer
    double metricBound = delta; // target metric bound
    double min_metric = -1;
    bool layerComplete;
    uint64_t min_candidate = 0;
    while(min_metric == -1) {
        // Growing forward tree
        for(unsigned int layer = 0; layer < n; ++layer) {
            // Adding eligible nodes from the buffer to the tree
            /*for(auto iter = nodeBuffer[layer].begin(); iter != nodeBuffer[layer].end(); ++iter)
            {
                if(iter->metric < (metricBound+delta))
                {
                    insertNode(*iter, fwdTree[layer]);
                    nodeBuffer[layer].erase(iter);
                }
            }*/
            // Getting layerMask for current layer
            layerMask = 0;
            layerComplete = false;
            for (unsigned int i = 0; i < k; ++i) {
                if (layer < ranges[2 * i] || layer >= ranges[2 * i + 1]) {
                    layerMask ^= (1 << i);
                }
            }
            for (auto iter = fwdTree[layer].begin(); iter != fwdTree[layer].end(); ++iter) {
                // If transition with input 0 is possible, add node to the tree
                if (layer < n &&
                    iter->path0 &&
                    !(iter->number & layerMask) &&
                    (iter->metric < metricBound)) {
                    temp.number = iter->number;
                    temp.metric = iter->metric + metric(x[layer], -1);
                    temp.path = iter->path;
                    temp.path0 = true;
                    temp.path1 = true;
                    // If resulting metric is too big, store it in buffer until metric bound is bigger
                    /*if(temp.metric > metricBound) {
                        nodeBuffer[layer+1].push_back(temp);
                    }
                    else {*/
                        insertNode(temp, fwdTree[layer+1]);
                    //}
                    iter->path0 = false;
                }
                // Same for input 1
                if (layer < n &&
                    iter->path1 &&
                    !((iter->number ^ h[layer]) & layerMask) &&
                    iter->metric < metricBound) {
                    temp.number = iter->number ^ h[layer];
                    temp.metric = iter->metric + metric(x[layer], 1);
                    temp.path = iter->path;
                    temp.path ^= (1 << layer);
                    temp.path0 = true;
                    temp.path1 = true;
                    /*if(temp.metric > metricBound) {
                        nodeBuffer[layer+1].push_back(temp);
                    }
                    else {*/
                        insertNode(temp, fwdTree[layer+1]);
                    //}
                    iter->path1 = false;
                }
                layerComplete = layerComplete || iter->path0 || iter->path1;
            }
            // If a tree layer does not have any nodes that can be continued, erase it
            if(!layerComplete)
            {
                fwdTree[layer].clear();
            }
        }
        // Growing backward tree
        for(int layer = n; layer > 0; --layer) {
            layerMask = 0;
            for (unsigned int i = 0; i < k; ++i) {
                if ((layer - 1) <= ranges[2 * i] || (layer - 1) > ranges[2 * i + 1]) {
                    layerMask ^= (1 << i);
                }
            }
            for (auto iter = bkwTree[layer].begin(); iter != bkwTree[layer].end(); ++iter) {
                // The difference between forward tree is in third condition
                if (layer > 0 &&
                    iter->path0 &&
                    !(iter->number & layerMask) &&
                    (iter->metric + metric(x[layer - 1], -1)) < metricBound) {
                    temp.number = iter->number;
                    temp.metric = iter->metric + metric(x[layer - 1], -1);
                    temp.path = iter->path;
                    temp.path0 = true;
                    temp.path1 = true;
                    insertNode(temp, bkwTree[layer-1]);
                    iter->path0 = false;
                }
                if (layer > 0 &&
                    iter->path1 &&
                    !((iter->number ^ h[layer - 1]) & layerMask) &&
                    (iter->metric + metric(x[layer - 1], 1)) < metricBound) {
                    temp.number = iter->number ^ h[layer - 1];
                    temp.metric = iter->metric + metric(x[layer - 1], 1);
                    temp.path = iter->path;
                    temp.path ^= (1 << (layer - 1));
                    temp.path0 = true;
                    temp.path1 = true;
                    insertNode(temp, bkwTree[layer-1]);
                    iter->path1 = false;
                }
            }
        }
        // Looking for matches in sorted lists
        for(unsigned int layer = 0; layer <=n; ++layer) {
            auto fwdIter = fwdTree[layer].begin();
            auto bkwIter = bkwTree[layer].begin();
            double tempMetric;
            NodeCompare nodecmpr;
            while (fwdIter != fwdTree[layer].end() && bkwIter != bkwTree[layer].end()) {
                if (!nodecmpr(*fwdIter, *bkwIter) && !nodecmpr(*bkwIter, *fwdIter)) {
                    tempMetric = fwdIter->metric + bkwIter->metric;
                    if (min_metric == -1 || tempMetric < min_metric) {
                        // Found a match, storing it
                        min_metric = tempMetric;
                        min_candidate = fwdIter->path + bkwIter->path;
                    }
                    ++fwdIter;
                    ++bkwIter;
                } else if (nodecmpr(*fwdIter, *bkwIter)) {
                    ++fwdIter;
                } else {
                    ++bkwIter;
                }
            }
        }
        if(min_metric > 2*metricBound) // found metric can be higher than the actual minimum, double-checking
        {
            metricBound = min_metric / 2;
            min_metric = -1;
        }
        else
        {
            metricBound += delta;
        }
    }
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = (min_candidate & (1 << i)) ? 1 : 0;
    }

    delete [] fwdTree;
    delete [] bkwTree;
    delete [] nodeBuffer;
    return min_metric;
}
