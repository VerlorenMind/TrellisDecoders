#ifndef BEAST_TRELLISDECODER_H
#define BEAST_TRELLISDECODER_H

#include <fstream>
#include <cstring>

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
    void init(unsigned int n, unsigned int k, uint64_t* checkmatrix);
public:
    unsigned int op_add, op_cmp;
    uint64_t* h;
    TrellisDecoder(unsigned int n, unsigned int k, std::ifstream& filename);
    TrellisDecoder(unsigned int n, unsigned int k, uint64_t* checkmatrix);
    virtual double decode(double* x, unsigned int* u, double delta);
    ~TrellisDecoder();
};
#endif //BEAST_TRELLISDECODER_H
