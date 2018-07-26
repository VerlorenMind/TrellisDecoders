#ifndef BEAST_BEASTDECODER_H
#define BEAST_BEASTDECODER_H

#include <fstream>

struct Node {
    unsigned int number;
    unsigned int layer;
    double metric;
    std::vector<unsigned int> path;
    bool path0;
    bool path1;
};

class BeastDecoder {
private:
    unsigned int n;
    unsigned int k;
    unsigned int* h;
    unsigned int* ranges;

public:
    BeastDecoder(unsigned int n, unsigned int k, const char* matrix_file);
    void decode(double* x, unsigned int* u, double delta);
    ~BeastDecoder();
};

unsigned int* readMatrix(const char* filename, unsigned int n, unsigned int k);
void minspan_form(unsigned int n, unsigned int k, unsigned int* a);
unsigned int* find_ranges(unsigned int n, unsigned int k, unsigned int* a);

#endif //BEAST_BEASTDECODER_H
