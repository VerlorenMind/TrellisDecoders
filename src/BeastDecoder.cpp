#include <iostream>
#include <list>
#include <cmath>
#include <cassert>

#include "BeastDecoder.h"

void BeastDecoder::init()
{
    fwdTree = new Node[trellisSize];
    fwdTreeSize = 0;
    fwdTreeBuffer = new Node[trellisSize];
    fwdTreeBufferSize = 0;
    bkwTree = new Node[trellisSize];
    bkwTreeSize = 0;
    bkwTreeBuffer = new Node[trellisSize];
    bkwTreeBufferSize = 0;
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, std::ifstream& filename) : TrellisDecoder(n, k, filename)
{
    init();
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, int **h) : TrellisDecoder(n, k, h)
{
    init();
}

BeastDecoder::~BeastDecoder()
{
    delete[] fwdTree;
    delete[] fwdTreeBuffer;
    delete[] bkwTree;
    delete[] bkwTreeBuffer;
}

InsertionStatus BeastDecoder::insertNode(const Node& node)
{
    InsertionStatus status;
    uint64_t trellisNum = 0;
    if(node.tree == FWD)
    {
        int ind = 0;
        for (unsigned int i = 0; i < k; ++i)
        {
            if ((node.layer - 1) >= ranges[2 * i] && (node.layer - 1) < ranges[2 * i + 1])
            {
                trellisNum ^= (((node.number >> i) & 1) << ind);
                ++ind;
            }
        }
    }
    else
    {
        int ind = 0;
        for (unsigned int i = 0; i < k; ++i)
        {
            if ((node.layer) > ranges[2 * i] && (node.layer) <= ranges[2 * i + 1])
            {
                trellisNum ^= (((node.number >> i) & 1) << ind);
                ++ind;
            }
        }
    }
    assert(trellisNum < trellisProfile[node.layer]);
    if(trellis[node.layer][trellisNum].tree == NIL) {
        trellis[node.layer][trellisNum] = node;
        status = INSERTED;
    }
    else if(trellis[node.layer][trellisNum].tree == node.tree)
    {
        if(trellis[node.layer][trellisNum].metric > node.metric)
        {
            trellis[node.layer][trellisNum] = node;
            status = REPLACED;
        }
        else
        {
            status = DISCARDED;
        }
        ++op_cmp;
    }
    else
    {
        double metric = trellis[node.layer][trellisNum].metric + node.metric;
        ++op_add;
        uint64_t candidate = node.path + trellis[node.layer][trellisNum].path;
        if(min_metric == -1)
        {
            min_metric = metric;
            min_candidate = candidate;
            status = ENDED;
        }
        else
        {
            if(metric < min_metric)
            {
                min_metric = metric;
                min_candidate = candidate;
                if(candidate == min_candidate)
                {
                    status = REPLACED;
                }
                else
                {
                    status = ENDED;
                }
            }
            else
            {
                if(candidate == min_candidate)
                {
                    status = DISCARDED;
                }
                else
                {
                    status = ENDED;
                }
            }
            ++op_cmp;
        }
    }
    return status;
}

double BeastDecoder::decode(double *x, int *u, double delta)
{
    op_add = 0;
    op_cmp = 0;

    for(unsigned int i=0; i<n; ++i)
    {
        alpha[i] = x[i] < 0 ? 0 : 1;
        ++op_cmp;
        beta[i] = fabs(x[i]);
    }

    // Emptying trees TODO: replace with memset
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
    while(true) {
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
                        ++op_cmp;
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
                        uint64_t tempnum = iter.number;
                        for(unsigned int i=0; i<k; ++i)
                        {
                            tempnum ^= (uint64_t(h[i][iter.layer]) << i);
                        }
                        if (tempnum & layerMask) {
                            iter.pathAvalaible[1] = false;
                        }

                        // Expanding tree
                        for (int k = 0; k<2; ++k) {
                            if (iter.pathAvalaible[k]) {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer+1;
                                temp.number = k ? tempnum : iter.number;
                                temp.metric = iter.metric+metric(k, iter.layer);
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
                                temp.pathAvalaible[0] = true;
                                temp.pathAvalaible[1] = true;
                                InsertionStatus status = insertNode(temp);
                                bool shouldBeContinued = temp.metric < fwdMetricBound;
                                ++op_cmp;
                                switch(status)
                                {
                                    case INSERTED:
                                    {
                                        if (shouldBeContinued)
                                        {
                                            assert(fwdTreeSize < trellisSize);
                                            fwdTree[fwdTreeSize] = temp;
                                            ++fwdTreeSize;
                                        } else
                                        {
                                            unsigned bufferIndex = 0;
                                            while (fwdTreeBuffer[bufferIndex].tree != NIL)
                                            {
                                                ++bufferIndex;
                                            }
                                            fwdTreeBuffer[bufferIndex] = temp;
                                            ++fwdTreeBufferSize;
                                        }
                                        ++fwdSize;
                                        break;
                                    }
                                    case REPLACED:
                                    {
                                        if (shouldBeContinued)
                                        {
                                            unsigned int newIndex = index, newCount = count;
                                            bool isInQueue = false;
                                            while(newIndex < trellisSize && newCount < fwdTreeSize)
                                            {
                                                Node newIter = fwdTree[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if(newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                        fwdTree[newIndex] = temp;
                                                        isInQueue = true;
                                                        break;
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                            if(!isInQueue)
                                            {
                                                assert(fwdTreeSize < trellisSize);
                                                fwdTree[fwdTreeSize] = temp;
                                                ++fwdTreeSize;
                                            }
                                        }
                                        unsigned int newIndex = 0;
                                        unsigned int newCount = 0;
                                        if(fwdTreeBufferSize > 0)
                                        {
                                            while (newIndex < trellisSize && newCount < fwdTreeBufferSize)
                                            {
                                                Node newIter = fwdTreeBuffer[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if (newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                        if (shouldBeContinued)
                                                        {
                                                            fwdTreeBuffer[newIndex].tree = NIL;
                                                            --fwdTreeBufferSize;
                                                        }
                                                        else
                                                        {
                                                            fwdTreeBuffer[newIndex] = temp;
                                                        }
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                        }
                                    }
                                    case DISCARDED:
                                    case ENDED: break;
                                }
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
                            bkwTreeBuffer[index].tree = NIL;
                            InsertionStatus status = insertNode(iter);
                            switch(status)
                            {
                                case INSERTED:
                                {
                                    bkwTree[bkwTreeSize] = iter;
                                    ++bkwSize;
                                    ++bkwTreeSize;
                                    break;
                                }
                                case REPLACED:
                                case DISCARDED:
                                case ENDED: break;
                            }
                        }
                        ++op_cmp;
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
                        uint64_t tempnum = iter.number;
                        for(unsigned int i=0; i<k; ++i)
                        {
                            tempnum ^= (uint64_t(h[i][iter.layer-1]) << i);
                        }
                        if (tempnum & layerMask) {
                            iter.pathAvalaible[1] = false;
                        }
                        for (int k = 0; k<2; ++k) {
                            double calcMetric = metric(k, iter.layer-1);
                            if (iter.layer>0 && iter.pathAvalaible[k])
                            {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer - 1;
                                temp.number = k ? tempnum : iter.number;
                                temp.metric = iter.metric + calcMetric;
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << (iter.layer - 1)) : iter.path;
                                temp.pathAvalaible[0] = true;
                                temp.pathAvalaible[1] = true;
                                bool shouldBeContinued = temp.metric < bkwMetricBound;
                                ++op_cmp;
                                if (shouldBeContinued)
                                {
                                    InsertionStatus status = insertNode(temp);
                                    switch (status)
                                    {
                                        case INSERTED:
                                        case ENDED:
                                        {
                                            assert(bkwTreeSize < trellisSize);
                                            bkwTree[bkwTreeSize] = temp;
                                            ++bkwTreeSize;
                                            ++bkwSize;
                                            break;
                                        }
                                        case REPLACED:
                                        {
                                            unsigned int newIndex = index, newCount = count;
                                            bool isInQueue = false;
                                            while (newIndex < trellisSize && newCount < bkwTreeSize)
                                            {
                                                Node newIter = bkwTree[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if (newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                        bkwTree[newIndex] = temp;
                                                        isInQueue = true;
                                                        break;
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                            if (!isInQueue)
                                            {
                                                assert(bkwTreeSize < trellisSize);
                                                bkwTree[bkwTreeSize] = temp;
                                                ++bkwTreeSize;
                                            }
                                            newIndex = 0;
                                            newCount = 0;
                                            if (bkwTreeBufferSize > 0)
                                            {
                                                while (newIndex < trellisSize && newCount < bkwTreeBufferSize)
                                                {
                                                    Node newIter = bkwTreeBuffer[newIndex];
                                                    if (newIter.tree != NIL)
                                                    {
                                                        ++newCount;
                                                        if (newIter.layer == temp.layer &&
                                                            newIter.number == temp.number)
                                                        {
                                                            bkwTreeBuffer[newIndex].tree = NIL;
                                                            --bkwTreeBufferSize;
                                                            break;
                                                        }
                                                    }
                                                    ++newIndex;
                                                }
                                                break;
                                            }
                                        }
                                        case DISCARDED:
                                            break;
                                    }
                                } else
                                {
                                    unsigned bufferIndex = 0;
                                    while (bkwTreeBuffer[bufferIndex].tree != NIL)
                                    {
                                        ++bufferIndex;
                                    }
                                    bkwTreeBuffer[bufferIndex] = temp;
                                    ++bkwTreeBufferSize;
                                }
                            }
                        }
                    }
                }
            }
        }
        if(min_metric != -1)
        {
            if(min_metric <= fwdMetricBound + bkwMetricBound)
            {
                ++op_cmp;
                break;
            }
            ++op_cmp;
        }
    }
    // Outputting the result
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = (min_candidate & (uint64_t(1) << i)) ? 1 : 0;
    }

    return min_metric;
}
