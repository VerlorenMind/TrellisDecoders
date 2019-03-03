#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include <cmath>
#include <BeastDecoder.h>

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
            if(a[i] & (uint64_t(1)<<j))
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
                for(unsigned int l=0; l<n; ++l)
                {
                    temp = a[l] & (uint64_t(1) << rows[0]);
                    a[l] ^= (((a[l] & (uint64_t(1) << i)) << (rows[0]-i)) ^ (a[l] & (uint64_t(1) << rows[0])));
                    a[l] ^= (temp >> (rows[0]-i)) ^ (a[l] & (uint64_t(1)<<i));
                }
                rows[0] = i;
            }

            ++fixed_rows;

            // Adding fixed row to the unfixed rest of the matrix
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[j] ^= (a[j] & (uint64_t(1)<<rows[0]))<<(rows[l]-rows[0]);
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
            if(a[i] & (uint64_t(1)<<j))
            {
                rows[num++] = (unsigned int) j;
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
                    a[j] ^= (a[j] & (uint64_t(1)<<rows[0]))>>(rows[0]-rows[l]);
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
                ranges[2*i+1] = (unsigned int) j;
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
    alpha = new int[n];
    beta = new double[n];
    h = readMatrix(filename, n, k);
    minspan_form(n, k, h);
    ranges = find_ranges(n, k, h);
    fwdTree = new std::set<Node, NodeCompare>[n+1];
    bkwTree = new std::set<Node, NodeCompare>[n+1];
    bkwTreeBuffer = new std::set<Node, MetricCompare>[n+1];
    maxLayerSize = k < n-k ? (uint64_t(2)<<k) : (uint64_t(2)<<(n-k));
}

BeastDecoder::~BeastDecoder()
{
    delete[] h;
    delete[] ranges;
    delete[] alpha;
    delete[] beta;
    delete[] fwdTree;
    delete[] bkwTree;
    delete[] bkwTreeBuffer;
}

inline double BeastDecoder::metric(int x, unsigned int pos)
{
    ++op_cmp;
    return (x == alpha[pos] ? 0 : beta[pos]);
}
template <class T>
void BeastDecoder::insertNode(const Node& node, std::set<Node, T>& tree)
{
    // Check if this node is already in the tree
    auto tempfind = tree.find(node);
    // If it isn't, add it
    ++op_cmp;
    if (tempfind == tree.end()) {
        tree.insert(node);
    } else {
        // If it is, try to minimize it's metric
        ++op_cmp;
        if (tempfind->metric > node.metric) {
            tempfind->metric = node.metric;
            tempfind->path = node.path;
        }
    }
}

double BeastDecoder::decode(double *x, unsigned int *u, double delta)
{
    op_add = 0;
    op_cmp = 0;
    op_mul = 0;
    op_bit = 0;
    for(unsigned int i=0; i<n; ++i)
    {
        alpha[i] = x[i] < 0 ? 0 : 1;
        beta[i] = fabs(x[i]);
        fwdTree[i].clear();
        bkwTree[i].clear();
        bkwTreeBuffer[i].clear();
    }
    fwdTree[n].clear();
    bkwTree[n].clear();
    bkwTreeBuffer[n].clear();
    // Initializing starting nodes
    Node temp{};
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.pathAvalaible[0] = true;
    temp.pathAvalaible[1] = true;
    fwdTree[0].insert(temp);
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.pathAvalaible[0] = true;
    temp.pathAvalaible[1] = true;
    bkwTree[n].insert(temp);
    uint64_t layerMask; // a mask that denotes zero bit-positions on a layer
    double fwdMetricBound = delta, bkwMetricBound = delta; // target metric bound
    double min_metric = -1;
    bool layerComplete;
    unsigned fwdSize = 1;
    unsigned bkwSize = 1;
    uint64_t min_candidate = 0;
    char metricChanged = 3;
    while(min_metric == -1) {
        // Growing forward tree
        if(metricChanged & 1) {
            for (unsigned int layer = 0; layer < n; ++layer) {
                // Getting layerMask for current layer
                layerMask = 0;
                layerComplete = false;
                for (unsigned int i = 0; i < k; ++i) {
                    if (layer < ranges[2 * i] || layer >= ranges[2 * i + 1]) {
                        layerMask ^= (uint64_t(1) << i);
                    }
                }
                for (auto iter = fwdTree[layer].begin(); iter != fwdTree[layer].end(); ++iter) {
                    if (iter->number & layerMask) {
                        iter->pathAvalaible[0] = false;
                    }
                    if ((iter->number ^ h[layer]) & layerMask) {
                        iter->pathAvalaible[1] = false;
                    }
                    if(iter->metric < fwdMetricBound) {
                        for (int k = 0; k < 2; ++k) {
                            if (layer < n &&
                                iter->pathAvalaible[k]) {
                                temp.number = k ? (iter->number ^ h[layer]) : iter->number;
                                temp.metric = iter->metric + metric(k, layer);
                                temp.path = k ? iter->path ^ (uint64_t(1) << layer) : iter->path;
                                temp.pathAvalaible[0] = true;
                                temp.pathAvalaible[1] = true;
                                insertNode<NodeCompare>(temp, fwdTree[layer + 1]);
                                ++fwdSize;
                                iter->pathAvalaible[k] = false;
                            }
                        }
                    }
                    layerComplete = layerComplete || iter->pathAvalaible[0] || iter->pathAvalaible[1];
                }
                // If a tree layer does not have any nodes that can be continued, erase it
                if (!layerComplete) {
                    fwdTree[layer].clear();
                }
                if (bkwTree[layer].size() > maxLayerSize) {
                    std::cerr << "FwdTree layer " << layer << " overflow: "
                              << fwdTree[layer].size() << " nodes when " << maxLayerSize << " should be\n"
                              << "For target metric: " << fwdMetricBound << std::endl;
                    return -1;
                }
            }
        }
        // Growing backward tree
        if(metricChanged & 2) {

            for (int layer = n; layer > 0; --layer) {
                layerMask = 0;
                temp.metric = bkwMetricBound;
                auto buffer_iter = bkwTreeBuffer[layer].upper_bound(temp);
                while(buffer_iter != bkwTreeBuffer[layer].end())
                {
                    insertNode<NodeCompare>(*buffer_iter, bkwTree[layer]);
                    bkwTreeBuffer[layer].erase(buffer_iter++);
                    ++bkwSize;
                }
                for (unsigned int i = 0; i < k; ++i) {
                    if ((layer - 1) <= ranges[2 * i] || (layer - 1) > ranges[2 * i + 1]) {
                        layerMask ^= (uint64_t(1) << i);
                    }
                }
                for (auto iter = bkwTree[layer].begin(); iter != bkwTree[layer].end(); ++iter) {
                    // The difference between forward tree is in third condition
                    if (iter->number & layerMask) {
                        iter->pathAvalaible[0] = false;
                    }
                    if ((iter->number ^ h[layer - 1]) & layerMask) {
                        iter->pathAvalaible[1] = false;
                    }
                    for (int k = 0; k < 2; ++k) {
                        double calcMetric = metric(k, (unsigned int) layer - 1);
                        if (layer > 0 && iter->pathAvalaible[k]) {
                            temp.number = k ? (iter->number ^ h[layer - 1]) : iter->number;
                            temp.metric = iter->metric + calcMetric;
                            temp.path = k ? iter->path ^ (uint64_t(1) << (layer - 1)) : iter->path;
                            temp.pathAvalaible[0] = true;
                            temp.pathAvalaible[1] = true;
                            if (temp.metric < bkwMetricBound) {
                                insertNode<NodeCompare>(temp, bkwTree[layer - 1]);
                                ++bkwSize;
                            }
                            else
                            {
                                insertNode<MetricCompare>(temp, bkwTreeBuffer[layer - 1]);
                            }
                            iter->pathAvalaible[k] = false;
                        }
                    }
                }
                if (bkwTree[layer].size() > maxLayerSize) {
                    std::cerr << "BkwTree layer " << layer << " overflow: "
                              << bkwTree[layer].size() << " nodes when " << maxLayerSize << " should be\n"
                              << "For target metric: " << bkwMetricBound << std::endl;
                    return -1;
                }
            }
        }

        // Looking for matches in sorted lists
        NodeCompare nodecmpr;
        for(unsigned int layer = 0; layer <=n; ++layer) {
            auto fwdIter = fwdTree[layer].begin();
            auto bkwIter = bkwTree[layer].begin();
            double tempMetric;
            // Iterating on both trees, incrementing iterators with lesser nodes
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
        // Found metric can be higher than the actual minimum, double-checking
        if(min_metric > fwdMetricBound+bkwMetricBound)
        {
            min_metric = -1;
        }
        if(bkwSize < fwdSize)
        {
            bkwMetricBound += delta;
            metricChanged = 2;
        }
        else
        {
            fwdMetricBound += delta;
            metricChanged = 1;
        }
    }
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = (min_candidate & (uint64_t(1) << i)) ? 1 : 0;
    }

    return min_metric;
}
