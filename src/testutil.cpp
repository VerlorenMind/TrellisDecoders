#include "testutil.h"

const double EPSILON = 0.0000001;

void generate_vector(unsigned int size, int* vector)
{
    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(0,1);
    for(unsigned int i=0; i<size; ++i)
    {
        vector[i] = distribution(generator);
    }
}

void test_decoder(unsigned int tests, double stn, double delta, const char* filename, std::string decoder)
{
    std::ifstream input(filename);
    unsigned int n, k;
    input >> n >> k;
    WARN("Test BCH(" << n << ", " << k << ") starts");
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    double overall = 0;

    double dev = sqrt((1 / (((1.0 * k) / n) * pow(10, stn / 10))) / 2);
    unsigned int seed = 123; // (unsigned int) std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::normal_distribution<double> rv(0.0, dev);

    int **g = readMatrix(input, n, k);
    int **h = readMatrix(input, n, n-k);
    TrellisDecoder *dec;
    if (decoder == "beast")
    {
        dec = new BeastDecoder(n, n-k, h);
    } else // if (decoder == "bsd")
    {
        dec = new BSDDecoder(n, n - k, h);
    }
    double* x = new double[n];
    int* y = new int[n];
    int* u = new int[k];
    int* ux = new int[n];
    double metric, calcWeight, noise, euclTrue, euclCalc;
    unsigned int temp;
    bool flag;
    std::string word;
    int fails = 0;
    int *synd = new int[n-k];
    for(unsigned int test = 0; test < tests; ++test)
    {
        generate_vector(k, u);
        INFO("Test #" << test);
        INFO("Informational word: " << array_to_sstream<int>(k, u).str());
        for(unsigned int i=0; i<n; ++i)
        {
            temp = 0;
            for(unsigned int j=0; j<k; ++j)
            {
                temp ^= g[j][i] & u[j];
            }
            ux[i] = (temp ? 1 : 0);
        }
        INFO("Coded word: " << array_to_sstream<int>(n, ux).str());
        for(unsigned int i=0; i<n; ++i)
        {
            noise = rv(gen);
            x[i] = (ux[i] ? 1. : -1.) + noise;
        }
        memset(synd, 0, (n-k)*sizeof(int));
        int syndsum = 0;
        for(unsigned int j=0; j<n-k; ++j)
        {
            for(unsigned int i=0; i<n; ++i)
            {
                synd[j] ^= ux[i] & h[j][i];
            }
        }
        for(unsigned int j=0; j<n-k; ++j)
        {
            syndsum += synd[j];
        }
        INFO("Coded word syndrome: " << array_to_sstream<int>(n-k, synd).str());
        REQUIRE(syndsum == 0);
        INFO("Coded word with noise: " << array_to_sstream<double>(n, x).str());
        metric = 0;
        for(unsigned int i=0; i<n; ++i)
        {
            if(x[i] < 0)
            {
                metric += (ux[i] == 0 ? 0 : fabs(x[i]));
            }
            else
            {
                metric += (ux[i] == 1 ? 0 : fabs(x[i]));
            }
        }
        euclTrue = 0;
        for(unsigned int i=0; i<n; ++i)
        {
            euclTrue += ((ux[i] ? 1. : -1.) - x[i])*((ux[i] ? 1. : -1.) - x[i]);
        }
        start = std::chrono::high_resolution_clock::now();

        calcWeight = dec->decode(x, y, delta);

        stop = std::chrono::high_resolution_clock::now();
        overall += std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        flag = true;
        for(unsigned int i=0; i<n; ++i)
        {
            if(y[i] != ux[i])
            {
                flag = false;
                break;
            }
        }
        uint64_t yint = 0;
        for(unsigned int i=0; i<n; ++i)
        {
            yint += y[i] << i;
        }
        memset(synd, 0, (n-k)*sizeof(int));
        syndsum = 0;
        for(unsigned int j=0; j<n-k; ++j)
        {
            for(unsigned int i=0; i<n; ++i)
            {
                synd[j] ^= y[i] & h[j][i];
            }
        }
        for(unsigned int j=0; j<n-k; ++j)
        {
            syndsum += synd[j];
        }
        euclCalc = 0;
        for(unsigned int i=0; i<n; ++i)
        {
            euclCalc += ((y[i] ? 1. : -1.) - x[i])*((y[i] ? 1. : -1.) - x[i]);
        }
        INFO("Decoded word: " << array_to_sstream<int>(n, y).str());
        INFO("True weight: " << metric);
        INFO("Calculated weight: " << calcWeight);
        INFO("Syndrome: " << array_to_sstream<int>(n-k, synd).str());
        INFO("True Euclidean metric: " << euclTrue);
        INFO("Resulted Euclidean metric: " << euclCalc);
        // INFO("Check matrix:\n"<<matrix_to_sstream(n-k, n, h).str());
        if(!((flag || calcWeight < metric || fabs(calcWeight - metric) < EPSILON) &&
            syndsum == 0 &&
            (euclCalc < euclTrue || fabs(euclCalc - euclTrue) < EPSILON)))
        {
            ++fails;
            CHECK(false);
        }
        if(!((test+1) % 250))
        {
            WARN("Time to decode "<< test+1<< " words: "<< overall<< "ms");
        }
    }

    delete[] x;
    delete[] y;
    delete[] u;
    for(unsigned int i=0; i<k; ++i)
    {
        delete[] g[i];
    }
    delete[] g;
    for(unsigned int i=0; i<n-k; ++i)
    {
        delete[] h[i];
    }
    delete[] h;
    delete[] ux;
    delete dec;
    delete[] synd;
    WARN("Overall time: "<< overall <<"ms");
    WARN("Failed tests: "<<fails);
}
