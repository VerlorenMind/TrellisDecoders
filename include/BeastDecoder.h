#ifndef BEAST_BEASTDECODER_H
#define BEAST_BEASTDECODER_H

#include <fstream>
#include <cstdint>
#include <set>
#include <queue>

enum Tree {
    NIL,
    FWD,
    BKW,
};

struct Node {
    Tree tree;
    unsigned layer;
    uint64_t number;
    mutable double metric;
    mutable uint64_t path;
    mutable bool pathAvalaible[2];
};

class NodeCompare {
public:
    unsigned int timesCalled = 0;
    bool operator()(const Node& lhs, const Node& rhs)
    {
        ++timesCalled;
        return lhs.number < rhs.number;
    }
};

class MetricCompare {
public:
    bool operator()(const Node& lhs, const Node& rhs)
    {
        return lhs.metric >= rhs.metric;
    }
};

class BeastDecoder {
private:
    unsigned int n;
    unsigned int k;
    unsigned int* ranges;
    int *alpha;
    double *beta;
    double min_metric = -1;
    uint64_t min_candidate = 0;
    uint64_t maxLayerSize;
    std::priority_queue<Node, std::vector<Node>, MetricCompare> fwdTree, bkwTree, bkwTreeBuffer;
    Node** trellis;
    unsigned* offsets;
    uint64_t* trellisProfile;
    // std::set<Node, NodeCompare> *fwdTree, *bkwTree;
    // std::set<Node, MetricCompare> *bkwTreeBuffer;

    void insertNode(const Node& node);
    inline double metric(int x, unsigned int pos);
public:
    unsigned int op_add, op_mul, op_cmp, op_bit;
    uint64_t* h;
    BeastDecoder(unsigned int n, unsigned int k, std::ifstream& filename);
    double decode(double* x, unsigned int* u, double delta);
    ~BeastDecoder();
};

uint64_t* readMatrix(std::ifstream& filename, unsigned int n, unsigned int k);
void minspan_form(unsigned int n, unsigned int k, uint64_t* a);
unsigned int* find_ranges(unsigned int n, unsigned int k, uint64_t* a);

#endif //BEAST_BEASTDECODER_H
