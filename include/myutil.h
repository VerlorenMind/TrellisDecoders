#ifndef BEAST_UTIL_H
#define BEAST_UTIL_H

#include <sstream>
#include <fstream>

template<class T>
std::stringstream array_to_sstream(unsigned int size, T* a)
{
    std::stringstream result;
    for(unsigned int i=0; i<size; ++i)
    {
        result << a[i] <<" ";
    }
    return result;
}

uint64_t* readMatrix(std::ifstream& input, unsigned int n, unsigned int k);

void minspan_form(unsigned int n, unsigned int k, uint64_t* a);

unsigned int* find_ranges(unsigned int n, unsigned int k, uint64_t* a);

#endif //BEAST_UTIL_H
