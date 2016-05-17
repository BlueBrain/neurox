#pragma once

#include <vector>

using namespace std;

/**
 * @brief The Mechanisms class
 * Stores unique mechanisms information and dependencies
 */
class Mechanism
{
  public:
    Mechanism() = delete;

    Mechanism( short int datasize, short int pdataSize,
        short int dependenciesCount, char pntMap,
        char isArtificial, int * dependencies);

    ~Mechanism();

    short int dataSize, pdataSize, dependenciesCount;
    char pntMap, isArtificial;

    //Dependencies data
    int * dependencies;

  private:
};
