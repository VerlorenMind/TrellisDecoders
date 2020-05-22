#include <iostream>
#include <list>
#include <cmath>
#include <cassert>

#include "BeastDecoder.h"

void BeastDecoder::init(double delta)
{
  id = DecoderID::BEAST;
  this->delta = delta;
  fwd_tree = new Node[trellis_size];
  fwd_tree_size = 0;
  fwd_tree_buffer = new Node[trellis_size];
  fwd_tree_buffer_size = 0;
  bkw_tree = new Node[trellis_size];
  bkw_tree_size = 0;
  bkw_tree_buffer = new Node[trellis_size];
  bkw_tree_buffer_size = 0;

    trellis = new Node *[n + 1];
    trellis[0] = new Node[1];
    for (unsigned i = 0; i < n; ++i)
    {
        trellis[i + 1] = new Node[trellis_profile[i + 1]];
        for (unsigned j = 0; j < trellis_profile[i]; ++j)
        {
            trellis[i][j].tree = NIL;
        }
    }
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, std::ifstream& filename, double delta) : TrellisDecoder(n, k, filename)
{
    init(delta);
}

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, int **h, double delta) : TrellisDecoder(n, k, h)
{
    init(delta);
}

BeastDecoder::~BeastDecoder()
{
    delete[] fwd_tree;
    delete[] fwd_tree_buffer;
    delete[] bkw_tree;
    delete[] bkw_tree_buffer;
    for(unsigned i=0; i<=n; ++i)
    {
        delete[] trellis[i];
    }
    delete[] trellis;
}

InsertionStatus BeastDecoder::insert_node(const Node& node)
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
    assert(trellisNum < trellis_profile[node.layer]);
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
                if(candidate == min_candidate)
                {
                    status = REPLACED;
                }
                else
                {
                    min_candidate = candidate;
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

double BeastDecoder::decode(double *x, int *u)
{
    op_add = 0;
    op_cmp = 0;

    for(unsigned int i=0; i<n; ++i)
    {
        alpha[i] = x[i] < 0 ? 0 : 1;
        ++op_cmp;
        beta[i] = fabs(x[i]);
    }

    // Emptying trees
    memset(fwd_tree, 0, trellis_size*sizeof(Node));
    memset(fwd_tree_buffer, 0, trellis_size*sizeof(Node));
    memset(bkw_tree, 0, trellis_size*sizeof(Node));
    memset(bkw_tree_buffer, 0, trellis_size*sizeof(Node));
    for(unsigned i=0; i<=n;++i)
    {
        for(unsigned j=0; j<trellis_profile[i]; ++j)
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
    temp.path_avalaible[0] = true;
    temp.path_avalaible[1] = true;
  fwd_tree[0] = temp;
  fwd_tree_size = 1;
    trellis[0][0] = temp;
    temp.tree = BKW;
    temp.layer = n;
  bkw_tree[0] = temp;
  bkw_tree_size = 1;
    trellis[n][0] = temp;
  fwd_tree_buffer_size = bkw_tree_buffer_size = 0;

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
            if(fwd_tree_buffer_size > 0) {
              fwd_tree_size = 0;
                while (index < trellis_size && count < fwd_tree_buffer_size) {
                    Node iter = fwd_tree_buffer[index];
                    if (iter.tree!=NIL) {
                        ++count;
                        if (iter.metric<fwdMetricBound) {
                          fwd_tree[fwd_tree_size] = iter;
                            ++fwd_tree_size;
                          fwd_tree_buffer[index].tree = NIL;
                        }
                        ++op_cmp;
                    }
                    ++index;
                }
              fwd_tree_buffer_size = fwd_tree_buffer_size - fwd_tree_size;
            }
            count = index = 0;
            while(index < trellis_size && count < fwd_tree_size) {
                Node iter = fwd_tree[index];
              fwd_tree[index].tree = NIL;
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
                            iter.path_avalaible[0] = false;
                        }
                        uint64_t tempnum = iter.number;
                        for(unsigned int i=0; i<k; ++i)
                        {
                            tempnum ^= (uint64_t(h[i][iter.layer]) << i);
                        }
                        if (tempnum & layerMask) {
                            iter.path_avalaible[1] = false;
                        }

                        // Expanding tree
                        for (int k = 0; k<2; ++k) {
                            if (iter.path_avalaible[k]) {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer+1;
                                temp.number = k ? tempnum : iter.number;
                                temp.metric = iter.metric+metric(k, iter.layer);
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
                                temp.path_avalaible[0] = true;
                                temp.path_avalaible[1] = true;
                                InsertionStatus status = insert_node(temp);
                                bool shouldBeContinued = temp.metric < fwdMetricBound;
                                ++op_cmp;
                                switch(status)
                                {
                                    case INSERTED:
                                    case ENDED:
                                    {
                                        if (shouldBeContinued)
                                        {
                                            assert(fwd_tree_size < trellis_size);
                                          fwd_tree[fwd_tree_size] = temp;
                                            ++fwd_tree_size;
                                        } else
                                        {
                                            unsigned bufferIndex = 0;
                                            while (fwd_tree_buffer[bufferIndex].tree != NIL)
                                            {
                                                ++bufferIndex;
                                            }
                                          fwd_tree_buffer[bufferIndex] = temp;
                                            ++fwd_tree_buffer_size;
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
                                            while(newIndex < trellis_size && newCount < fwd_tree_size)
                                            {
                                                Node newIter = fwd_tree[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if(newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                      fwd_tree[newIndex] = temp;
                                                        isInQueue = true;
                                                        break;
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                            if(!isInQueue)
                                            {
                                                assert(fwd_tree_size < trellis_size);
                                              fwd_tree[fwd_tree_size] = temp;
                                                ++fwd_tree_size;
                                            }
                                        }
                                        unsigned int newIndex = 0;
                                        unsigned int newCount = 0;
                                        if(fwd_tree_buffer_size > 0)
                                        {
                                            while (newIndex < trellis_size && newCount < fwd_tree_buffer_size)
                                            {
                                                Node newIter = fwd_tree_buffer[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if (newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                        if (shouldBeContinued)
                                                        {
                                                          fwd_tree_buffer[newIndex].tree = NIL;
                                                            --fwd_tree_buffer_size;
                                                        }
                                                        else
                                                        {
                                                          fwd_tree_buffer[newIndex] = temp;
                                                        }
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                        }
                                    }
                                    case DISCARDED: break;
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
            if(bkw_tree_buffer_size > 0) {
              bkw_tree_size = 0;
                while (index < trellis_size && count < bkw_tree_buffer_size) {
                    Node iter = bkw_tree_buffer[index];
                    if (iter.tree!=NIL) {
                        ++count;
                        if (iter.metric<bkwMetricBound) {
                          bkw_tree_buffer[index].tree = NIL;
                            InsertionStatus status = insert_node(iter);
                            switch(status)
                            {
                                case INSERTED:
                                case ENDED:
                                {
                                  bkw_tree[bkw_tree_size] = iter;
                                    ++bkwSize;
                                    ++bkw_tree_size;
                                    break;
                                }
                                case REPLACED:
                                case DISCARDED: break;
                            }
                        }
                        ++op_cmp;
                    }
                    ++index;
                }
              bkw_tree_buffer_size = bkw_tree_buffer_size - bkw_tree_size;
            }
            count = index = 0;
            while(index < trellis_size && count < bkw_tree_size) {
                Node iter = bkw_tree[index];
              bkw_tree[index].tree = NIL;
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
                            iter.path_avalaible[0] = false;
                        }
                        uint64_t tempnum = iter.number;
                        for(unsigned int i=0; i<k; ++i)
                        {
                            tempnum ^= (uint64_t(h[i][iter.layer-1]) << i);
                        }
                        if (tempnum & layerMask) {
                            iter.path_avalaible[1] = false;
                        }
                        for (int k = 0; k<2; ++k) {
                            double calcMetric = metric(k, iter.layer-1);
                            if (iter.layer>0 && iter.path_avalaible[k])
                            {
                                temp.tree = iter.tree;
                                temp.layer = iter.layer - 1;
                                temp.number = k ? tempnum : iter.number;
                                temp.metric = iter.metric + calcMetric;
                                ++op_add;
                                temp.path = k ? iter.path ^ (uint64_t(1) << (iter.layer - 1)) : iter.path;
                                temp.path_avalaible[0] = true;
                                temp.path_avalaible[1] = true;
                                bool shouldBeContinued = temp.metric < bkwMetricBound;
                                ++op_cmp;
                                if (shouldBeContinued)
                                {
                                    InsertionStatus status = insert_node(temp);
                                    switch (status)
                                    {
                                        case INSERTED:
                                        case ENDED:
                                        {
                                            assert(bkw_tree_size < trellis_size);
                                          bkw_tree[bkw_tree_size] = temp;
                                            ++bkw_tree_size;
                                            ++bkwSize;
                                            break;
                                        }
                                        case REPLACED:
                                        {
                                            unsigned int newIndex = index, newCount = count;
                                            bool isInQueue = false;
                                            while (newIndex < trellis_size && newCount < bkw_tree_size)
                                            {
                                                Node newIter = bkw_tree[newIndex];
                                                if (newIter.tree != NIL)
                                                {
                                                    ++newCount;
                                                    if (newIter.layer == temp.layer && newIter.number == temp.number)
                                                    {
                                                      bkw_tree[newIndex] = temp;
                                                        isInQueue = true;
                                                        break;
                                                    }
                                                }
                                                ++newIndex;
                                            }
                                            if (!isInQueue)
                                            {
                                                assert(bkw_tree_size < trellis_size);
                                              bkw_tree[bkw_tree_size] = temp;
                                                ++bkw_tree_size;
                                            }
                                            newIndex = 0;
                                            newCount = 0;
                                            if (bkw_tree_buffer_size > 0)
                                            {
                                                while (newIndex < trellis_size && newCount < bkw_tree_buffer_size)
                                                {
                                                    Node newIter = bkw_tree_buffer[newIndex];
                                                    if (newIter.tree != NIL)
                                                    {
                                                        ++newCount;
                                                        if (newIter.layer == temp.layer &&
                                                            newIter.number == temp.number)
                                                        {
                                                          bkw_tree_buffer[newIndex].tree = NIL;
                                                            --bkw_tree_buffer_size;
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
                                    while (bkw_tree_buffer[bufferIndex].tree != NIL)
                                    {
                                        ++bufferIndex;
                                    }
                                  bkw_tree_buffer[bufferIndex] = temp;
                                    ++bkw_tree_buffer_size;
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
void BeastDecoder::set_delta(double d) {
  delta = d;
}
double BeastDecoder::get_delta() {
  return delta;
}
