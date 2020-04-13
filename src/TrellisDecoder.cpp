#include "TrellisDecoder.h"

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, std::ifstream &filename)
{
    int **checkmatrix = readMatrix(filename, n, k);
    init(n, k, checkmatrix);
    for(unsigned int i=0; i<k; ++i)
    {
        delete[] checkmatrix[i];
    }
    delete[] checkmatrix;
}

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, int **checkmatrix)
{
    init(n, k, checkmatrix);
}

TrellisDecoder::~TrellisDecoder()
{
    for(unsigned int i=0; i<k; ++i)
    {
       delete[] h[i];
    }
    delete[] h;
    delete[] ranges;
    delete[] alpha;
    delete[] beta;
    // delete[] offsets;
    for(unsigned i=0; i<=n; ++i)
    {
        delete[] trellis[i];
    }
    delete[] trellis;
    delete[] trellisProfile;
}

double TrellisDecoder::metric(int x, unsigned int pos)
{
    ++op_cmp;
    return (x == alpha[pos] ? 0 : beta[pos]);
}

double TrellisDecoder::decode(double *x, int *u, double delta)
{
    return 0;
}

void TrellisDecoder::init(unsigned int n, unsigned int k, int **checkmatrix)
{
// Transforming matrix to minspan form and finding active ranges of it's rows
    this->n = n;
    this->k = k;
    alpha = new int[n];
    beta = new double[n];
    h = new int*[k];
    for(unsigned int i=0; i<k; ++i)
    {
        h[i] = new int[n];
        memcpy(h[i], checkmatrix[i], n*sizeof(int));
    }
    std::string tempmatr = matrix_to_sstream(k, n, h).str();
    minspan_form(n, k, h);
    tempmatr = matrix_to_sstream(k, n, h).str();
    ranges = find_ranges(n, k, h);
    trellis = new Node *[n + 1];
    trellisProfile = new uint64_t[n + 1];
    unsigned pow;
    trellis[0] = new Node[1];
    trellisProfile[0] = 1;
    trellisSize = 1;
    for (unsigned i = 0; i < n; ++i)
    {
        pow = 0;
        for (unsigned int j = 0; j < k; ++j)
        {
            if (i >= ranges[2 * j] && i < ranges[2 * j + 1])
            {
                ++pow;
            }
        }
        trellisProfile[i + 1] = uint64_t(1) << pow;
        trellisSize += trellisProfile[i + 1];
        trellis[i + 1] = new Node[trellisProfile[i + 1]];
        for (unsigned j = 0; j < trellisProfile[i]; ++j)
        {
            trellis[i][j].tree = NIL;
        }
    }
    maxLayerSize = k < n-k ? (uint64_t(2)<<k) : (uint64_t(2)<<(n-k));
}
