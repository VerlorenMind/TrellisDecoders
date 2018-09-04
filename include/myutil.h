#ifndef BEAST_UTIL_H
#define BEAST_UTIL_H

#include "BeastDecoder.h"
#include <sstream>
#include <iostream>

template<class T>
void print_matrix(unsigned int n, unsigned int k, T** m)
{
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            std::cout<<m[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}

template<class T>
void print_vector(unsigned int n, T* v)
{
    for(unsigned int j=0; j<n; ++j)
    {
        std::cout<<v[j]<<" ";
    }
    std::cout<<"\n";
}

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

#endif //BEAST_UTIL_H
