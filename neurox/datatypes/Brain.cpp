#include "neurox/Neurox.h"
#include <cstring>

#include "coreneuron/nrnoc/multicore.h" //nrn_threads

Brain * brain = nullptr; //global variable (defined in Neurox.h)

Brain::~Brain()
{
    delete [] mechanisms;
}

Brain::Brain(const int neuronsCount,
             const hpx_t neuronsAddr, const Mechanism * mechanisms,
             const size_t mechanismsCount, int * mechDependencies)
    : neuronsCount(neuronsCount), neuronsAddr(neuronsAddr), mechanismsCount(mechanismsCount)
{
    //add mechanisms information
    int offset=0;
    this->mechanisms = new Mechanism[mechanismsCount];
    for (int m=0; m<mechanismsCount; m++)
    {
        this->mechanisms[m]=Mechanism(m, mechanisms[m].dataSize, mechanisms[m].pdataSize,
                                      mechanisms[m].dependenciesCount, mechanisms[m].pntMap,
                                      mechanisms[m].isArtificial, &mechDependencies[offset],
                                      mechanisms[m].isIon, mechanisms[m].conci, mechanisms[m].conco, mechanisms[m].charge //for ions
                                      );
        offset += mechanisms[m].dependenciesCount;
    }

    //copy BA functions
    for (int bat=0; bat<BEFORE_AFTER_SIZE; bat++)
    {
        int idMechFunction = bat+4;
        for (int nt=0; nt<nrn_nthread; nt++)
        {
            for (NrnThreadBAList *tbl = nt->tbl[bat]; tbl; tbl = tbl->next)
            {
                auto f = tbl->bam->f;
                int type = tbl->type;
                mechanisms[type].functions[idMechFunction]=f;
            }
        }
    }
}
