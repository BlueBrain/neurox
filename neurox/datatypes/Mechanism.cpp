#include "neurox/neurox.h"
#include <cstring>

using namespace std;

Mechanism::Mechanis( short int datasize, short int pdataSize,
                     short int dependenciesCount, char pntMap,
                     char isArtificial, int &* dependencies):
    this->dataSize(datasize), this->pdatasize(pdatasize),
    this->dependenciesCount(dependenciesCount),
    this->pntMap(pntMap), this->isArtificial(isArtificial)
{
    dependencies = new int[dependenciesCount];
    std::memcpy(this->dependencies, dependencies, dependenciesCount*sizeof(int));
};

Mechanism::~Mechanism(){
    delete [] dependencies;
};
