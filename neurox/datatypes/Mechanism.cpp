#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Mechanism::Mechanism(short int dataSize, short int pdataSize,
                     short int dependenciesCount, char pntMap,
                     char isArtificial, int * dependencies):
    dataSize(dataSize), pdataSize(pdataSize), dependenciesCount(dependenciesCount),
    pntMap(pntMap), isArtificial(isArtificial)
{
    dependencies = new int[dependenciesCount];
    std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
};

Mechanism::~Mechanism(){
    delete [] dependencies;
};
