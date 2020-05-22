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
    delete[] trellis_profile;
}

double TrellisDecoder::metric(int x, unsigned int pos)
{
    ++op_cmp;
    return (x == alpha[pos] ? 0 : beta[pos]);
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
  minspanForm(n, k, h);
    ranges = findRanges(n, k, h);
  trellis_profile = new uint64_t[n + 1];
    unsigned pow;
  trellis_profile[0] = 1;
  trellis_size = 1;
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
      trellis_profile[i + 1] = uint64_t(1) << pow;
      trellis_size += trellis_profile[i + 1];
    }
}
