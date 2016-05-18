#pragma once

#include <vector>
//#include "coreneuron/nrnoc/membfunc.h"

using namespace std;

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Mechanism
{
  public:
    Mechanism(){};
    ~Mechanism();

    Mechanism( const short int datasize, const short int pdataSize,
        const short int dependenciesCount, const char pntMap,
        const char isArtificial, const int * dependencies);

    short int dataSize, pdataSize, dependenciesCount;
    char pntMap, isArtificial;

    int * dependencies; ///> Id of dependencies mechanisms

  private:
};
