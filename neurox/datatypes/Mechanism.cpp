#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Mechanism::Mechanism(const short int dataSize, const short int pdataSize,
                     const short int dependenciesCount, const char pntMap,
                     const char isArtificial, const int * dependencies):
    dataSize(dataSize), pdataSize(pdataSize), dependenciesCount(dependenciesCount),
    pntMap(pntMap), isArtificial(isArtificial)
{
    dependencies = new int[dependenciesCount];
    std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
};

Mechanism::~Mechanism(){
    delete [] dependencies;
};
