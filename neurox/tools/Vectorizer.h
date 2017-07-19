#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>
#include <deque>

using namespace std;

namespace neurox {

namespace tools {

/**
 * @brief The Vectorizer class
 * converts code from AoS to SoA structure
 */
class Vectorizer
{
  public:
    Vectorizer()=delete;
    ~Vectorizer()=delete;

    ///get size of datastructure with alignment-based padding
    static size_t SizeOf(size_t);

    ///converts a branch from AoS to SoA
    static void ConvertToSOA(Branch * b);

    //C++11 does not support memory-aligned new[]/delete, this is a work around
    template<typename T>
    static T* New(size_t count)
    {
        return (T*)coreneuron::ecalloc_align(SizeOf(count), NEUROX_MEM_ALIGNMENT, sizeof(T));
    }

    template<typename T>
    static void Delete(T * ptr)
    {
        if (ptr==nullptr) return;
        //delete[] (ptr);
        free(ptr);
        ptr=nullptr;
    }

  private:
};

}; //Tools
}; //NeuroX

