#ifndef BEAST_UTIL_H
#define BEAST_UTIL_H

#include "BeastDecoder.h"
#include <sstream>
#include <iostream>

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
