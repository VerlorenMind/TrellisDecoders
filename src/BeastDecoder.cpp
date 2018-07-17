#define TESTING

#include <algorithm>
#include <iostream>
#include <vector>
#include <stack>
#include <list>
#include "myutil.h"
#include "BeastDecoder.h"

unsigned int* readMatrix(const char* filename, unsigned int n, unsigned int k)
{
    // Reading matrix from a file and storing it's columns as uints
    std::ifstream input(filename);
    unsigned int* matrix = new unsigned int[n];
    for(unsigned int i=0; i<n; ++i)
    {
        matrix[i] = 0;
    }
    unsigned int bit;
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            input >> bit;
            matrix[j] ^= (bit << i);
        }
    }
    return matrix;
}

void minspan_form(unsigned int n, unsigned int k, unsigned int* a) {
    unsigned int num; // number of rows with 1s in i'th column
    unsigned int fixed_rows = 0; // number of rows on top of the matrix that were transformed
    unsigned int* rows = new unsigned int[k]; // rows with 1s in i'th column
    // Left side
    for(unsigned int i=0; i<k; ++i)
    {
        num = 0;
        // Finding all rows with 1s in i'th column
        for(unsigned int j=fixed_rows; j<k; ++j)
        {
            if(a[i] & (1<<j))
            {
                rows[num++] = j;
            }
        }
        if(num == 0)
        {
            continue; // should never happen with non-zero determinant TODO: handle this exception
        }
        else
        {
            // If the first row after fixed rows does not contain 1, swap it with the first row that does
            if(rows[0] != fixed_rows)
            {
                unsigned int temp;
                for(unsigned int l=0; l<n; ++l)
                {
                    // Swapping bits
                    temp = a[l] & (1 << rows[0]);
                    a[l] ^= (((a[l] & (1 << i)) << (rows[0]-i)) ^ (a[l] & (1 << rows[0])));
                    a[l] ^= (temp >> (rows[0]-i)) ^ (a[l] & (1<<i));
                }
                rows[0] = i;
            }

            ++fixed_rows;

            // Adding fixed row to the unfixed rest of the matrix
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[j] ^= (a[j] & (1<<rows[0]))<<(rows[l]-rows[0]);
                }
            }
        }
    }
    // Right side
    // Same stuff as above, but with different indices and without swapping rows
    fixed_rows = 0;
    for(unsigned int i=n-1; i>=n-k; --i)
    {
        num = 0;
        for(int j=k-fixed_rows-1; j>=0; --j)
        {
            if(a[i] & (1<<j))
            {
                rows[num++] = j;
            }
        }
        if(num == 0)
        {
            continue;
        }
        else
        {
            ++fixed_rows;
            for(unsigned int l=1; l<num; ++l)
            {
                for(unsigned int j=0; j<n; ++j)
                {
                    a[j] ^= (a[j] & (1<<rows[0]))>>(rows[0]-rows[l]);
                }
            }
        }
    }
    delete [] rows;
}

unsigned int* find_ranges(unsigned int n, unsigned int k, unsigned int* a)
{
    // Finding active ranges for each row (range between the first non-zero element and the last one)
    unsigned int* ranges = new unsigned int[2*k];
    for(unsigned int i=0; i<k; ++i)
    {
        for(unsigned int j=0; j<n; ++j)
        {
            if(a[j] & (1<<i))
            {
                ranges[2*i] = j;
                break;
            }
        }
        for(unsigned int j=n-1; j>=0; --j) {
            if (a[j] & (1<<i)) {
                ranges[2*i+1] = j;
                break;
            }
        }
    }
    return ranges;
}

/*Node** construct_trellis(unsigned int n, unsigned int k, bool** h)
{
    std::list<unsigned int> active_symb;
    Node **nodelayers = new Node*[n+1];
    nodelayers[0] = new Node[1];
    int *ranges = new int[2*k];
    find_ranges(n, k, h, ranges);
    nodelayers[0][0].number = 0;
    nodelayers[0][0].metric = 0;
    nodelayers[0][0].layer = 0;
    nodelayers[0][0].pathfwd0 = nullptr;
    nodelayers[0][0].pathfwd1 = nullptr;
    nodelayers[0][0].pathbkw0 = nullptr;
    nodelayers[0][0].pathbkw1 = nullptr;
    unsigned int tempnum;
    unsigned int curnum, xornum;
    unsigned int prevsize = 1;
    unsigned int newsize;
    unsigned int ind;
    unsigned int act_num;
    active_symb.push_back(0);
    for(unsigned int j=0; j<n; ++j)
    {
        newsize = 1<<(active_symb.size());
        nodelayers[j+1] = new Node[newsize];
        tempnum = 0;
        for(unsigned int i=0; i<newsize; ++i)
        {
            act_num = 0;
            ind = 0;
            for(auto active = active_symb.begin(); active != active_symb.end(); ++active)
            {
                act_num ^= ((tempnum & (1 << ind)) >> ind) << (*active);
                ++ind;
            }
            nodelayers[j+1][i].number = act_num;
            nodelayers[j+1][i].layer = j + 1;
            nodelayers[j+1][i].metric = -1;
            nodelayers[j+1][i].pathfwd0 = nullptr;
            nodelayers[j+1][i].pathfwd1 = nullptr;
            nodelayers[j+1][i].pathbkw0 = nullptr;
            nodelayers[j+1][i].pathbkw1 = nullptr;
            tempnum += 1;
        }
        for(unsigned int i=0; i<prevsize; ++i)
        {
            curnum = nodelayers[j][i].number;
            xornum = curnum;
            for(unsigned int l=0; l<k; ++l)
            {
                xornum ^= ((h[l][j] ? 1 : 0) << l);
            }
            for(unsigned int l=0; l<newsize; ++l)
            {
                if(nodelayers[j+1][l].number == curnum)
                {
                    nodelayers[j][i].pathfwd0 = nodelayers[j+1] + l;
                    nodelayers[j+1][l].pathbkw0 = nodelayers[j] + i;
                }
                if(nodelayers[j+1][l].number == xornum)
                {
                    nodelayers[j][i].pathfwd1 = nodelayers[j+1] + l;
                    nodelayers[j+1][l].pathbkw1 = nodelayers[j] + i;
                }
            }
        }
        active_symb.clear();
        for(unsigned int i=0; i<k; ++i)
        {
            if((j+1) >= ranges[2*i] && (j+1) < ranges[2*i+1])
            {
                active_symb.push_back(i);
            }
        }
        prevsize = newsize;
    }
    nodelayers[n][0].metric = 0;
    delete [] ranges;
    return nodelayers;
} */

BeastDecoder::BeastDecoder(unsigned int n, unsigned int k, const char* matrix_file)
{
    this->n = n;
    this->k = k;
    h = readMatrix(matrix_file, n, k);
    minspan_form(n, k, h);
    ranges = find_ranges(n, k, h);
}

BeastDecoder::~BeastDecoder()
{
    delete[] h;
    delete[] ranges;
}

inline double metric(double x, bool y)
{
    return y ? (x-1)*(x-1) : x*x;
}

void BeastDecoder::decode(double *x, unsigned int *u, double delta)
{

}

/*void printTree(Node start)
{
    std::deque<Node> outputDeque;
    outputDeque.push_back(start);
    bool flag;
    Node temp;
    while(!outputDeque.empty())
    {
        temp = outputDeque.front();
        outputDeque.pop_front();
        std::cout<<"Layer: "<<temp.layer<<"\n";
        std::cout<<"Number: "<<temp.number<<"\n";
        if(temp.pathfwd0 != nullptr) {
            std::cout << "0 -> " << (*temp.pathfwd0).number << "\n";
            flag = true;
            for(auto elem = outputDeque.begin(); elem != outputDeque.end(); ++elem)
            {
                if(elem->layer == (*temp.pathfwd0).layer && elem->number == (*temp.pathfwd0).number)
                {
                    flag = false;
                    break;
                }
            }
            if(flag) {
                outputDeque.push_back((*temp.pathfwd0));
            }
        }
        if(temp.pathfwd1 != nullptr) {
            std::cout << "1 -> " << (*temp.pathfwd1).number << "\n";
            flag = true;
            for(auto elem = outputDeque.begin(); elem != outputDeque.end(); ++elem)
            {
                if(elem->layer == (*temp.pathfwd1).layer && elem->number == (*temp.pathfwd1).number)
                {
                    flag = false;
                    break;
                }
            }
            if(flag) {
                outputDeque.push_back((*temp.pathfwd1));
            }
        }
        std::cout<<"---\n";
    }
}*/