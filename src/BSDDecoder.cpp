#include <cmath>
#include <cstring>
#include "BSDDecoder.h"

void BSDDecoder::bubbleInsert(Node* array, unsigned start, unsigned end, Node elem)
{
    if(start == end)
    {
        array[start] = elem;
    }
    unsigned i;
    for(i=start; i<end; ++i)
    {
        if(elem.metric < array[i].metric)
        {
            for(unsigned j=i; j<end+1; ++j)
            {
                std::swap(elem, array[j]);
            }
            break;
        }
        ++op_cmp;
    }
    if(i == end)
    {
        array[i] = elem;
    }
}

BSDDecoder::BSDDecoder(unsigned int n, unsigned int k, std::ifstream& filename) : TrellisDecoder(n, k, filename)
{
    fwdTree = new Node[2*trellisSize];
    bkwTree = new Node[2*trellisSize];
}


BSDDecoder::BSDDecoder(unsigned int n, unsigned int k, int **h) : TrellisDecoder(n, k, h)
{
    fwdTree = new Node[2*trellisSize];
    bkwTree = new Node[2*trellisSize];
}

BSDDecoder::~BSDDecoder()
{
    delete[] fwdTree;
    delete[] bkwTree;
}

double BSDDecoder::decode(double *x, int *u, double delta)
{
    op_add = 0;
    op_cmp = 0;

    for(unsigned int i=0; i<n; ++i)
    {
        alpha[i] = x[i] < 0 ? 0 : 1;
        ++op_cmp;
        beta[i] = fabs(x[i]);
    }
    memset(fwdTree, 0, 2*trellisSize*sizeof(Node));
    memset(bkwTree, 0, 2*trellisSize*sizeof(Node));

    min_metric = -1;
    uint64_t min_candidate = 0;
    Node temp{};
    temp.tree = FWD;
    temp.layer = 0;
    temp.number = 0;
    temp.metric = 0;
    temp.path = 0;
    temp.pathAvalaible[0] = true;
    temp.pathAvalaible[1] = true;
    fwdTree[0] = temp;
    temp.tree = BKW;
    temp.layer = n;
    bkwTree[0] = temp;

    int fwdLayerReached = 0;
    int bkwLayerReached = n;
    int fwdTreeStart = 0;
    int fwdTreeEnd = 1;
    int bkwTreeStart = 0;
    int bkwTreeEnd = 1;

    uint64_t layerMask;
    while(min_metric == -1)
    {
        Node iter;
        if(fwdTreeEnd - fwdTreeStart > 0)
        {
            iter = fwdTree[fwdTreeStart];
            fwdTree[fwdTreeStart].tree = NIL;
            if(iter.layer < n)
            {
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
                for (unsigned int j = 0; j < 2; ++j) {
                    // If path is not available, pathAvailable <- false and skipping
                    if (iter.pathAvalaible[j]) {
                        temp.tree = FWD;
                        temp.layer = iter.layer + 1;
                        fwdLayerReached = fwdLayerReached < temp.layer ? temp.layer : fwdLayerReached;
                        temp.number = j ? tempnum : iter.number;
                        double m = metric(j, iter.layer);
                        temp.path = j ? iter.path ^ (uint64_t(1) << iter.layer) : iter.path;
                        temp.pathAvalaible[0] = true;
                        temp.pathAvalaible[1] = true;
                        if (m != 0) {
                            temp.metric = iter.metric + m;
                            ++op_add;
                            bubbleInsert(fwdTree, fwdTreeStart + 1, fwdTreeEnd, temp);
                            ++fwdTreeEnd;
                        } else {
                            temp.metric = iter.metric;
                            fwdTree[fwdTreeStart] = temp;
                        }
                        ++op_cmp;
                    }
                }
                if (fwdTree[fwdTreeStart].tree == NIL && fwdTreeStart < fwdTreeEnd) {
                    ++fwdTreeStart;
                }
            }
        }
        if (bkwLayerReached - fwdLayerReached <= 0) {
            for (unsigned l = bkwTreeStart; l < bkwTreeEnd; ++l) {
                for (unsigned j = fwdTreeStart; j < fwdTreeEnd; ++j) {
                    if (fwdTree[j].layer == bkwTree[l].layer &&
                        fwdTree[j].number == bkwTree[l].number) {
                        double m = fwdTree[j].metric + bkwTree[l].metric;
                        ++op_add;
                        if (min_metric == -1) {
                            min_metric = m;
                            min_candidate = bkwTree[l].path + fwdTree[j].path;
                        }
                        else
                        {
                            if (m < min_metric)
                            {
                                min_metric = m;
                                min_candidate = bkwTree[l].path + fwdTree[j].path;
                            }
                            ++op_cmp;
                        }
                    }
                }
            }
        }
        if(bkwTreeEnd - bkwTreeStart > 0)
        {
            iter = bkwTree[bkwTreeStart];
            bkwTree[bkwTreeStart].tree = NIL;
            if(iter.layer > 0)
            {
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
                for (unsigned int j = 0; j < 2; ++j) {
                    // If path is not available, pathAvailable <- false and skipping
                    if (iter.pathAvalaible[j]) {
                        temp.tree = BKW;
                        temp.layer = iter.layer - 1;
                        bkwLayerReached = bkwLayerReached > temp.layer ? temp.layer : bkwLayerReached;
                        temp.number = j ? tempnum: iter.number;
                        double m = metric(j, iter.layer-1);
                        temp.path = j ? iter.path ^ (uint64_t(1) << (iter.layer - 1)) : iter.path;
                        temp.pathAvalaible[0] = true;
                        temp.pathAvalaible[1] = true;
                        if (m != 0) {
                            temp.metric = iter.metric + m;
                            ++op_add;
                            bubbleInsert(bkwTree, bkwTreeStart + 1, bkwTreeEnd, temp);
                            ++bkwTreeEnd;
                        } else {
                            temp.metric = iter.metric;
                            bkwTree[bkwTreeStart] = temp;
                        }
                        ++op_cmp;
                    }
                }
                if (bkwTree[bkwTreeStart].tree == NIL && bkwTreeStart < bkwTreeEnd) {
                    ++bkwTreeStart;
                }
            }
        }
        if (bkwLayerReached - fwdLayerReached <= 0) {
            for (unsigned l = fwdTreeStart; l < fwdTreeEnd; ++l) {
                for (unsigned j = bkwTreeStart; j < bkwTreeEnd; ++j) {
                    if (bkwTree[j].layer == fwdTree[l].layer &&
                        fwdTree[l].number == bkwTree[j].number) {
                        double m = fwdTree[l].metric + bkwTree[j].metric;
                        ++op_add;
                        if (min_metric == -1) {
                            min_metric = m;
                            min_candidate = bkwTree[j].path + fwdTree[l].path;
                        }
                        else
                        {
                            if (m < min_metric)
                            {
                                min_metric = m;
                                min_candidate = bkwTree[j].path + fwdTree[l].path;
                            }
                            ++op_cmp;
                        }
                    }
                }
            }
        }
    }
    for(unsigned int i=0; i<n; ++i)
    {
        u[i] = (min_candidate & (uint64_t(1) << i)) ? 1 : 0;
    }
    return min_metric;
}

