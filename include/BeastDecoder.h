#ifndef BEAST_BEASTDECODER_H
#define BEAST_BEASTDECODER_H

#include <fstream>
#include <cstdint>

struct Node {
    uint64_t number;
    unsigned int layer;
    double metric;
    uint64_t path;
    bool path0;
    bool path1;
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
