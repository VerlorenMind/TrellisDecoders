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
            // If the first row after fixed rows does not contain 1, swap it with the first row that does TODO: rewrite this part
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
    trellis = new Node*[n+1];
    trellisProfile = new uint64_t[n+1];
    offsets = new unsigned[n+1];
    bool offsetFlag;
    unsigned pow;
    trellis[0] = new Node[1];
    trellisProfile[0] = 1;
    offsets[0] = 0;
    for(unsigned i=0; i<n; ++i)
    {
        pow = 0;
        offsets[i+1] = 0;
        offsetFlag = true;
        for (unsigned int j = 0; j < k; ++j) {
            if (i >= ranges[2 * j] && i < ranges[2 * j + 1]) {
                ++pow;
                offsetFlag = false;
            }
            else if(offsetFlag)
            {
                ++offsets[i+1];
            }
        }
        trellisProfile[i+1] = uint64_t(1)<<pow;
        trellis[i+1] = new Node[trellisProfile[i+1]];
        for(unsigned j=0; j<trellisProfile[i]; ++j)
        {
            trellis[i][j].tree = NIL;
        }
    }
    // fwdTree = new std::set<Node, NodeCompare>[n+1];
    // bkwTree = new std::set<Node, NodeCompare>[n+1];
    // bkwTreeBuffer = new std::set<Node, MetricCompare>[n+1];
    maxLayerSize = k < n-k ? (uint64_t(2)<<k) : (uint64_t(2)<<(n-k));
}

BeastDecoder::~BeastDecoder()
{
    delete[] h;
    delete[] ranges;
    delete[] alpha;
    delete[] beta;
    delete[] offsets;
    for(unsigned i=0; i<=n; ++i)
    {
        delete[] trellis[i];
    }
    delete[] trellis;
    // delete[] fwdTree;
    // delete[] bkwTree;
    // delete[] bkwTreeBuffer;
}

inline double BeastDecoder::metric(int x, unsigned int pos)
{
    ++op_cmp;
    return (x == alpha[pos] ? 0 : beta[pos]);
}

void BeastDecoder::insertNode(const Node& node)
{
    uint64_t trellisNum = node.number >> offsets[node.layer];
    if(trellis[node.layer][trellisNum].tree == NIL) {
        trellis[node.layer][trellisNum] = node;
    }
    else if(trellis[node.layer][trellisNum].tree == node.tree)
    {
        if(trellis[node.layer][trellisNum].metric > node.metric)
        {
            trellis[node.layer][trellisNum] = node;
        }
        ++op_cmp;
    }
    else
    {
        double metric = trellis[node.layer][trellisNum].metric + node.metric;
        ++op_add;
        if(min_metric == -1 || metric < min_metric)
        {
            min_metric = metric;
            min_candidate = node.path + trellis[node.layer][trellisNum].path;
        }
        ++op_cmp;
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
    }

    // Emptying trees
    while(!fwdTree.empty())
    {
        fwdTree.pop();
    }
    while(!bkwTree.empty())
    {
        bkwTree.pop();
    }
    while(!bkwTreeBuffer.empty())
    {
        bkwTreeBuffer.pop();
    }
    for(unsigned i=0; i<=n;++i)
    {
        for(unsigned j=0; j<trellisProfile[i]; ++j)
        {
            trellis[i][j].tree = NIL;
        }
    }

    // Initializing starting nodes
    Node temp{};
    temp.tree = FWD;
    temp.layer = 0;
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.pathAvalaible[0] = true;
    temp.pathAvalaible[1] = true;
    fwdTree.push(temp);
    trellis[0][0] = temp;
    temp.tree = BKW;
    temp.layer = n;
    bkwTree.push(temp);
    trellis[n][0] = temp;

    uint64_t layerMask; // a mask that denotes zero bit-positions on a layer
    double fwdMetricBound = 0, bkwMetricBound = 0; // target metric bound
    unsigned fwdSize = 1;
    unsigned bkwSize = 1;
    min_metric = -1;
    min_candidate = 0;
    char metricChanged;
    while(min_metric == -1 || min_metric > fwdMetricBound+bkwMetricBound) {
        if(bkwSize < fwdSize)
        {
            bkwMetricBound += delta;
            ++op_add;
            metricChanged = 2;
        }
        else if(bkwSize > fwdSize)
        {
            fwdMetricBound += delta;
            ++op_add;
            metricChanged = 1;
        }
        else
        {
            fwdMetricBound += delta;
            bkwMetricBound += delta;
            op_add += 2;
            metricChanged = 3;
        }
        op_cmp += 2;
        ++op_add;
        // Growing forward tree
        if(metricChanged & 1) {
            Node iter = fwdTree.top();
            ++op_cmp;
            while (iter.metric < n && iter.metric < fwdMetricBound) {
                ++op_cmp;
                fwdTree.pop();

                // Getting layerMask for current layer TODO move these calculations to constructor
                layerMask = 0;
                for (unsigned int i = 0; i < k; ++i) {
                    if (iter.layer < ranges[2 * i] || iter.layer >= ranges[2 * i + 1]) {
                        layerMask ^= (uint64_t(1) << i);
                    }
                }

                // Checking avalaible continuations for a node
                if (iter.number & layerMask) {
                    iter.pathAvalaible[0] = false;
                }
                if ((iter.number ^ h[iter.layer]) & layerMask) {
                    iter.pathAvalaible[1] = false;
                }

                // Expanding tree
                for (int k = 0; k < 2; ++k) {
                    if (iter.layer < n &&
                        iter.pathAvalaible[k]) {
                        temp.tree = iter.tree;
                        temp.layer = iter.layer + 1;
                        temp.number = k ? (iter.number ^ h[iter.layer]) : iter.number;
                        temp.metric = iter.metric + metric(k, iter.layer);
                        ++op_add;
                        temp.path = k ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
                        temp.pathAvalaible[0] = true;
                        temp.pathAvalaible[1] = true;
                        insertNode(temp);
                        fwdTree.push(temp);
                        ++fwdSize;
                    }
                }
                iter = fwdTree.top();
            }
        }
        // Growing backward tree
        if(metricChanged & 2)
        {
            if(!bkwTreeBuffer.empty()) {
                Node buffer_iter = bkwTreeBuffer.top();
                while (buffer_iter.metric < bkwMetricBound) {
                    ++op_cmp;
                    bkwTreeBuffer.pop();
                    bkwTree.push(buffer_iter);
                    insertNode(buffer_iter);
                    ++bkwSize;
                    buffer_iter = bkwTreeBuffer.top();
                }
                ++op_cmp;
            }

            Node iter = bkwTree.top();
            ++op_cmp;
            while (iter.layer > 0 && iter.metric < bkwMetricBound) {
                ++op_cmp;

                if(bkwTree.empty())
                {
                    break;
                }
                bkwTree.pop();
                // std::cerr<<"Size after pop: "<<bkwTree.size()<<std::endl;

                // TODO: Move these computations to constructor
                layerMask = 0;
                for (unsigned int i = 0; i < k; ++i) {
                    if ((iter.layer - 1) <= ranges[2 * i] || (iter.layer - 1) > ranges[2 * i + 1]) {
                        layerMask ^= (uint64_t(1) << i);
                    }
                }

                if (iter.number & layerMask) {
                    iter.pathAvalaible[0] = false;
                }
                if ((iter.number ^ h[iter.layer - 1]) & layerMask) {
                    iter.pathAvalaible[1] = false;
                }
                for (int k = 0; k < 2; ++k) {
                    double calcMetric = metric(k, iter.layer - 1);
                    if (iter.layer > 0 && iter.pathAvalaible[k]) {
                        temp.tree = iter.tree;
                        temp.layer = iter.layer - 1;
                        temp.number = k ? (iter.number ^ h[iter.layer - 1]) : iter.number;
                        temp.metric = iter.metric + calcMetric;
                        ++op_add;
                        temp.path = k ? iter.path ^ (uint64_t(1) << (iter.layer - 1)) : iter.path;
                        temp.pathAvalaible[0] = true;
                        temp.pathAvalaible[1] = true;
                        if (temp.metric < bkwMetricBound) {
                            insertNode(temp);
                            bkwTree.push(temp);
                            // std::cerr<<"Size after push: "<<bkwTree.size()<<std::endl;
                            ++bkwSize;
                        }
                        else
                        {
                            bkwTreeBuffer.push(temp);
                        }
                        ++op_cmp;
                    }
                }
                if(bkwTree.empty())
                {
                    break;
                }
                iter = bkwTree.top();
            }
        }
    }
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = (min_candidate & (uint64_t(1) << i)) ? 1 : 0;
    }

    return min_metric;
}
