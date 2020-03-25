#ifndef BEAST_TRELLISDECODER_H
#define BEAST_TRELLISDECODER_H

#include <fstream>

#include "myutil.h"

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

class TrellisDecoder {
protected:
    unsigned int n;
    unsigned int k;
    unsigned int* ranges;
    int *alpha;
    double *beta;
    Node** trellis;
    unsigned* offsets;
    unsigned trellisSize;
    uint64_t* trellisProfile;
    uint64_t maxLayerSize;
    double metric(int x, unsigned int pos);
public:
    unsigned int op_add, op_mul, op_cmp, op_bit;
    uint64_t* h;
    TrellisDecoder(unsigned int n, unsigned int k, std::ifstream& filename);
    virtual double decode(double* x, unsigned int* u, double delta);
    ~TrellisDecoder();
};
#endif //BEAST_TRELLISDECODER_H
