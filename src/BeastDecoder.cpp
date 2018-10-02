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
            // If the first row after fixed rows does not contain 1, sum it with the first row that does
            if(rows[0] != fixed_rows)
            {
                for(unsigned int l=0; l<n; ++l)
                {
                    // a[l] ^= (a[l] & (uint64_t(1) << rows[0])) >> (rows[0] - fixed_rows);
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
    alpha = new int[n];
    beta = new double[n];
    h = readMatrix(filename, n, k);
    minspan_form(n, k, h);
    ranges = find_ranges(n, k, h);
    fwdTree = new std::set<Node, NodeCompare>[n+1];
    bkwTree = new std::set<Node, NodeCompare>[n+1];
}

BeastDecoder::~BeastDecoder()
{
    delete[] h;
    delete[] ranges;
    delete[] alpha;
    delete[] beta;
    delete[] fwdTree;
    delete[] bkwTree;
}

inline double BeastDecoder::metric(int x, unsigned int pos)
{
    ++op_cmp;
    return (x == alpha[pos] ? 0 : beta[pos]);
}

void BeastDecoder::insertNode(Node& node, std::set<Node, NodeCompare>& tree)
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
        ++op_cmp;
        ++op_add;
        alpha[i] = x[i] < 0 ? 0 : 1;
        ++op_cmp;
        beta[i] = fabs(x[i]);
        fwdTree[i].clear();
        bkwTree[i].clear();
    }
    fwdTree[n].clear();
    bkwTree[n].clear();
    //std::list<Node>* nodeBuffer = new std::list<Node>[n+1];
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
    unsigned int layerMask; // a mask that denotes zero bit-positions on a layer
    double metricBound = delta; // target metric bound
    double min_metric = -1;
    bool layerComplete;
    uint64_t min_candidate = 0;
    while(min_metric == -1) {
        ++op_cmp;
        // Growing forward tree
        std::cout<<"Target metric: "<<delta<<"\n"<<"\tForward tree:\n";

        for(unsigned int layer = 0; layer < n; ++layer) {
            ++op_cmp;
            ++op_add;
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
                ++op_cmp;
                ++op_add;
                op_cmp += 2;
                if (layer < ranges[2 * i] || layer >= ranges[2 * i + 1]) {
                    op_bit += 2;
                    layerMask ^= (1 << i);
                }
            }
            for (auto iter = fwdTree[layer].begin(); iter != fwdTree[layer].end(); ++iter) {
                ++op_cmp;
                ++op_add;
                ++op_bit;
                if(iter->number & layerMask)
                {
                    iter->pathAvalaible[0] = false;
                }
                op_cmp += 2;
                if (layer < n &&
                    iter->pathAvalaible[0] &&
                    (iter->metric < metricBound)) {
                    temp.number = iter->number;
                    temp.metric = iter->metric + metric(0, layer);
                    ++op_add;
                    temp.path = iter->path;
                    temp.pathAvalaible[0] = true;
                    temp.pathAvalaible[1] = true;
                    // If resulting metric is too big, store it in buffer until metric bound is bigger
                    //if(temp.metric > metricBound) {
                    //    nodeBuffer[layer+1].push_back(temp);
                    //}
                    //else {
                    insertNode(temp, fwdTree[layer + 1]);
                    //}
                    iter->pathAvalaible[0] = false;
                }
                op_bit += 2;
                if((iter->number ^ h[layer]) & layerMask)
                {
                    iter->pathAvalaible[1] = false;
                }
                op_cmp += 2;
                if (layer < n &&
                    iter->pathAvalaible[1] &&
                    (iter->metric < metricBound)) {
                    temp.number = (iter->number ^ h[layer]);
                    ++op_bit;
                    temp.metric = iter->metric + metric(1, layer);
                    ++op_add;
                    temp.path = iter->path;
                    temp.path ^= (1 << layer);
                    op_bit += 2;
                    temp.pathAvalaible[0] = true;
                    temp.pathAvalaible[1] = true;
                    // If resulting metric is too big, store it in buffer until metric bound is bigger
                    //if(temp.metric > metricBound) {
                    //    nodeBuffer[layer+1].push_back(temp);
                    //}
                    //else {
                    insertNode(temp, fwdTree[layer + 1]);
                    //}
                    iter->pathAvalaible[1] = false;
                }
                op_bit += 2;
                layerComplete = layerComplete || iter->pathAvalaible[0] || iter->pathAvalaible[1];
            }
            // If a tree layer does not have any nodes that can be continued, erase it
            ++op_bit;
            if(!layerComplete)
            {
                fwdTree[layer].clear();
            }
            std::cout<<"\t\tLayer "<<layer<<" nodes: "<<fwdTree[layer].size()<<"\n";
        }
        // Growing backward tree
        std::cout<<"\tBackward tree:\n";
        for(int layer = n; layer > 0; --layer) {
            ++op_cmp;
            ++op_add;
            layerMask = 0;
            for (unsigned int i = 0; i < k; ++i) {
                ++op_cmp;
                ++op_add;
                op_cmp += 2;
                op_mul += 2;
                op_add += 3;
                if ((layer - 1) <= ranges[2 * i] || (layer - 1) > ranges[2 * i + 1]) {
                    op_bit += 2;
                    layerMask ^= (1 << i);
                }
            }
            for (auto iter = bkwTree[layer].begin(); iter != bkwTree[layer].end(); ++iter) {
                ++op_cmp;
                ++op_add;
                // The difference between forward tree is in third condition
                op_cmp += 2;
                op_bit += 2;
                op_add += 2;
                if (layer > 0 &&
                    iter->pathAvalaible[0] &&
                    !(iter->number & layerMask) &&
                    (iter->metric + metric(0, layer-1)) < metricBound) {
                    temp.number = iter->number;
                    temp.metric = iter->metric + metric(0, layer-1);
                    op_add += 2;
                    temp.path = iter->path;
                    temp.pathAvalaible[0] = true;
                    temp.pathAvalaible[1] = true;
                    insertNode(temp, bkwTree[layer-1]);
                    ++op_add;
                    iter->pathAvalaible[0] = false;
                }
                op_cmp += 2;
                op_bit += 3;
                op_add += 1;
                if (layer > 0 &&
                    iter->pathAvalaible[1] &&
                    !((iter->number ^ h[layer - 1]) & layerMask) &&
                    (iter->metric + metric(1, layer-1)) < metricBound) {
                    temp.number = iter->number ^ h[layer - 1];
                    ++op_bit;
                    temp.metric = iter->metric + metric(1, layer-1);
                    ++op_add;
                    temp.path = iter->path;
                    temp.path ^= (1 << (layer - 1));
                    op_bit += 2;
                    ++op_add;
                    temp.pathAvalaible[0] = true;
                    temp.pathAvalaible[1] = true;
                    insertNode(temp, bkwTree[layer-1]);
                    iter->pathAvalaible[1] = false;
                }
            }
            std::cout<<"\t\tLayer "<<layer<<" nodes: "<<bkwTree[layer].size()<<"\n";
        }
        // Looking for matches in sorted lists
        NodeCompare nodecmpr;
        for(unsigned int layer = 0; layer <=n; ++layer) {
            ++op_cmp;
            ++op_add;
            auto fwdIter = fwdTree[layer].begin();
            auto bkwIter = bkwTree[layer].begin();
            double tempMetric;
            while (fwdIter != fwdTree[layer].end() && bkwIter != bkwTree[layer].end()) {
                op_bit += 3;
                if (!nodecmpr(*fwdIter, *bkwIter) && !nodecmpr(*bkwIter, *fwdIter)) {
                    ++op_add;
                    tempMetric = fwdIter->metric + bkwIter->metric;
                    op_cmp += 2;
                    if (min_metric == -1 || tempMetric < min_metric) {
                        // Found a match, storing it
                        min_metric = tempMetric;
                        ++op_add;
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
            op_cmp += nodecmpr.timesCalled;
        }
        ++op_cmp;
        ++op_mul;
        if(min_metric > 2*metricBound) // found metric can be higher than the actual minimum, double-checking
        {
            ++op_mul;
            metricBound = min_metric / 2;
            min_metric = -1;
        }
        else
        {
            ++op_add;
            metricBound += delta;
        }
    }
    op_cmp += fwdTree->key_comp().timesCalled + bkwTree->key_comp().timesCalled;
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        ++op_cmp;
        op_bit += 2;
        u[i] = (min_candidate & (1 << i)) ? 1 : 0;
    }

    //delete [] nodeBuffer;
    return min_metric;
}
