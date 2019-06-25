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
    trellisSize = 1;
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
        trellisSize += trellisProfile[i+1];
        trellis[i+1] = new Node[trellisProfile[i+1]];
        for(unsigned j=0; j<trellisProfile[i]; ++j)
        {
            trellis[i][j].tree = NIL;
        }
    }
    fwdTree = new Node[trellisSize];
    fwdTreeBuffer = new Node[trellisSize];
    bkwTree = new Node[trellisSize];
    bkwTreeBuffer = new Node[trellisSize];
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
    delete[] fwdTree;
    delete[] fwdTreeBuffer;
    delete[] bkwTree;
    delete[] bkwTreeBuffer;
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
    for(unsigned i=0; i<trellisSize; ++i)
    {
        fwdTree[i].tree = NIL;
        fwdTreeBuffer[i].tree = NIL;
        bkwTree[i].tree = NIL;
        bkwTreeBuffer[i].tree = NIL;
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
    fwdTree[0] = temp;
    fwdTreeSize = 1;
    trellis[0][0] = temp;
    temp.tree = BKW;
    temp.layer = n;
    bkwTree[0] = temp;
    bkwTreeSize = 1;
    trellis[n][0] = temp;
    fwdTreeBufferSize = bkwTreeBufferSize = 0;

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
            op_cmp += 1;
            bkwMetricBound += delta;
            ++op_add;
            metricChanged = 2;
        }
        else if(bkwSize > fwdSize)
        {
            op_cmp += 2;
            fwdMetricBound += delta;
            ++op_add;
            metricChanged = 1;
        }
        else
        {
            op_cmp += 2;
            fwdMetricBound += delta;
            bkwMetricBound += delta;
            op_add += 2;
            metricChanged = 3;
        }
        // Growing forward tree
        if(metricChanged & 1) {
            unsigned count = 0, index = 0;
            if(fwdTreeBufferSize > 0) {
                fwdTreeSize = 0;
                while (index < trellisSize && count < fwdTreeBufferSize) {
                    Node iter = fwdTreeBuffer[index];
                    if (iter.tree!=NIL) {
                        ++count;
                        if (iter.metric<fwdMetricBound) {
                            fwdTree[fwdTreeSize] = iter;
                            ++fwdTreeSize;
                            fwdTreeBuffer[index].tree = NIL;
                        }
                    }
                    ++index;
                }
                fwdTreeBufferSize = fwdTreeBufferSize - fwdTreeSize;
            }
            count = index = 0;
            while(index < trellisSize && count < fwdTreeSize) {
                Node iter = fwdTree[index];
                fwdTree[index].tree = NIL;
                ++index;
                if(iter.tree != NIL) {
                    ++count;
                    if(iter.layer < n) {
                        // Getting layerMask for current layer TODO move these calculations to constructor
                        layerMask = 0;
                        for (unsigned int i = 0; i<k; ++i) {
                            if (iter.layer<ranges[2*i] || iter.layer>=ranges[2*i+1]) {
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
                        for (int k = 0; k<2; ++k) {
                            if (iter.pathAvalaible[k]) {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer+1;
                                temp.number = k ? (iter.number ^ h[iter.layer]) : iter.number;
                                temp.metric = iter.metric+metric(k, iter.layer);
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
                                temp.pathAvalaible[0] = true;
                                temp.pathAvalaible[1] = true;
                                insertNode(temp);
                                if(temp.metric < fwdMetricBound) {
                                    fwdTree[fwdTreeSize] = temp;
                                    ++fwdTreeSize;
                                }
                                else {
                                    unsigned bufferIndex = 0;
                                    while(fwdTreeBuffer[bufferIndex].tree != NIL) {
                                        ++bufferIndex;
                                    }
                                    fwdTreeBuffer[bufferIndex] = temp;
                                    ++fwdTreeBufferSize;
                                }
                                ++fwdSize;
                            }
                        }
                    }
                }
            }
        }
        // Growing backward tree
        if(metricChanged & 2)
        {
            unsigned count = 0, index = 0;
            if(bkwTreeBufferSize > 0) {
                bkwTreeSize = 0;
                while (index < trellisSize && count < bkwTreeBufferSize) {
                    Node iter = bkwTreeBuffer[index];
                    if (iter.tree!=NIL) {
                        ++count;
                        if (iter.metric<bkwMetricBound) {
                            bkwTree[bkwTreeSize] = iter;
                            ++bkwTreeSize;
                            bkwTreeBuffer[index].tree = NIL;
                            insertNode(iter);
                        }
                    }
                    ++index;
                }
                bkwTreeBufferSize = bkwTreeBufferSize - bkwTreeSize;
            }
            count = index = 0;
            while(index < trellisSize && count < bkwTreeSize) {
                Node iter = bkwTree[index];
                bkwTree[index].tree = NIL;
                ++index;
                if(iter.tree != NIL) {
                    ++count;
                    if(iter.layer > 0) {
                        // TODO: Move these computations to constructor
                        layerMask = 0;
                        for (unsigned int i = 0; i<k; ++i) {
                            if ((iter.layer-1)<=ranges[2*i] || (iter.layer-1)>ranges[2*i+1]) {
                                layerMask ^= (uint64_t(1) << i);
                            }
                        }

                        if (iter.number & layerMask) {
                            iter.pathAvalaible[0] = false;
                        }
                        if ((iter.number ^ h[iter.layer-1]) & layerMask) {
                            iter.pathAvalaible[1] = false;
                        }
                        for (int k = 0; k<2; ++k) {
                            double calcMetric = metric(k, iter.layer-1);
                            if (iter.layer>0 && iter.pathAvalaible[k]) {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer-1;
                                temp.number = k ? (iter.number ^ h[iter.layer-1]) : iter.number;
                                temp.metric = iter.metric+calcMetric;
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << (iter.layer-1)) : iter.path;
                                temp.pathAvalaible[0] = true;
                                temp.pathAvalaible[1] = true;
                                if (temp.metric<bkwMetricBound) {
                                    insertNode(temp);
                                    bkwTree[bkwTreeSize] = temp;
                                    ++bkwTreeSize;
                                    ++bkwSize;
                                }
                                else {
                                    unsigned bufferIndex = 0;
                                    while (bkwTreeBuffer[bufferIndex].tree!=NIL) {
                                        ++bufferIndex;
                                    }
                                    bkwTreeBuffer[bufferIndex] = temp;
                                    ++bkwTreeBufferSize;
                                }
                                ++op_cmp;
                            }
                        }
                    }
                }
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
