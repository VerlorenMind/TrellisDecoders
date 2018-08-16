#ifndef BEAST_BEASTDECODER_H
#define BEAST_BEASTDECODER_H

#include <fstream>
#include <cstdint>

struct Node {
    uint64_t number;
    unsigned int layer;
    mutable double metric;
    mutable uint64_t path;
    mutable bool path0;
    mutable bool path1;
};

class NodeCompare {
public:
    bool operator()(const Node& lhs, const Node& rhs)
    {
        if(lhs.layer != rhs.layer) {
            return lhs.layer < rhs.layer;
        }
        else {
            return lhs.number < rhs.number;
        }
    }
};

class BeastDecoder {
private:
    unsigned int n;
    unsigned int k;
    uint64_t* h;
    unsigned int* ranges;

public:
    BeastDecoder(unsigned int n, unsigned int k, std::ifstream& filename);
    double decode(double* x, unsigned int* u, double delta);
    ~BeastDecoder();
};

uint64_t* readMatrix(std::ifstream& filename, unsigned int n, unsigned int k);
void minspan_form(unsigned int n, unsigned int k, uint64_t* a);
unsigned int* find_ranges(unsigned int n, unsigned int k, uint64_t* a);

#endif //BEAST_BEASTDECODER_H
