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

    ///converts a branch from AoS to SoA
    static void convertToSOA(Branch * b);

    //C++11 does not support memory-aligned new[]/delete, this is a work around
    template<typename T>
    static T* new_(size_t count)
    {
        void* ptr=nullptr;
        int err = posix_memalign(&ptr, NEUROX_MEM_ALIGNMENT, sizeof(T)*count);
        assert(err==0);
        //std::memset(ptr,0, sizeof(T)*count); return (T*) ptr;
        return new (ptr) T[count]();
    }

    template<typename T>
    static void delete_(T * ptr)
    {
        if (ptr==nullptr) return;
        delete[] (ptr); //free(ptr);
        ptr=nullptr;
    }

    template<typename T>
    static T* enlarge_(T* ptr, size_t count, size_t new_count)
    {
        assert(new_count>=count && new_count>0);
        T* new_ptr = new_<T>(new_count);
        assert(new_ptr && ptr);
        memcpy(new_ptr, ptr, sizeof(T)*count);
        return new_ptr;
    }

    ///get size of datastructure with alignment-based padding
    static size_t sizeof_(size_t);

  private:
};

}; //Tools
}; //NeuroX

