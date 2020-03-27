#ifndef BEAST_BSDDECODER_H
#define BEAST_BSDDECODER_H

#include "TrellisDecoder.h"

class BSDDecoder : public TrellisDecoder
{
protected:
    Node* fwdTree;
    Node* bkwTree;

    double min_metric;
    void bubbleInsert(Node* array, unsigned start, unsigned end, Node elem);
public:
    BSDDecoder(unsigned int n, unsigned int k, std::ifstream& filename);
    BSDDecoder(unsigned int n, unsigned int k, uint64_t *h); 
    ~BSDDecoder();
    double decode(double *x, unsigned int *u, double delta);
};
#endif //BEAST_BSDDECODER_H
