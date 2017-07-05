#pragma once

#include "neurox/neurox.h"

#include <map>
#include <vector>
#include <memory>
#include <deque>

using namespace std;

namespace neurox {

namespace Tools {

/**
 * @brief The Vectorizer class
 * converts code from AoS to SoA structure
 */
class Vectorizer
{
  public:
    Vectorizer()=delete;
    ~Vectorizer()=delete;

    static void vectorize(Branch * b);
};

}; //Tools
}; //NeuroX

