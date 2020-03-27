#include "TrellisDecoder.h"

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, std::ifstream &filename)
{
    uint64_t *checkmatrix = readMatrix(filename, n, k);
    init(n, k, checkmatrix);
    delete[] checkmatrix;
}

TrellisDecoder::TrellisDecoder(unsigned int n, unsigned int k, uint64_t* checkmatrix)
{
    init(n, k, checkmatrix);
}

TrellisDecoder::~TrellisDecoder()
{
    delete[] h;
    delete[] ranges;
    delete[] alpha;
    delete[] beta;
    delete[] offsets;
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

double TrellisDecoder::decode(double *x, unsigned int *u, double delta)
{
    return 0;
}

void TrellisDecoder::init(unsigned int n, unsigned int k, uint64_t *checkmatrix)
{
// Transforming matrix to minspan form and finding active ranges of it's rows
    this->n = n;
    this->k = k;
    alpha = new int[n];
    beta = new double[n];
    h = new uint64_t[n];
    memcpy(h, checkmatrix, n*sizeof(uint64_t));
    minspan_form(n, k, h);
    ranges = find_ranges(n, k, h);
    trellis = new Node *[n + 1];
    trellisProfile = new uint64_t[n + 1];
    offsets = new unsigned[n + 1];
    bool offsetFlag;
    unsigned pow;
    trellis[0] = new Node[1];
    trellisProfile[0] = 1;
    offsets[0] = 0;
    trellisSize = 1;
    for (unsigned i = 0; i < n; ++i)
    {
        pow = 0;
        offsets[i + 1] = 0;
        offsetFlag = true;
        for (unsigned int j = 0; j < k; ++j)
        {
            if (i >= ranges[2 * j] && i < ranges[2 * j + 1])
            {
                ++pow;
                offsetFlag = false;
            } else if (offsetFlag)
            {
                ++offsets[i + 1];
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
