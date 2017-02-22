#pragma once

#include "neurox/neurox.h"

#include <memory>
#include <vector>
#include <utility>

using namespace std;
using namespace neurox;
using namespace neurox::Input::Coreneuron;

class neurox::NetConX;

Compartment::~Compartment(){};

Compartment::Compartment(offset_t id, floble_t a, floble_t b, floble_t d,
                         floble_t v, floble_t rhs, floble_t area, offset_t p):
    id(id), a(a), b(b), d(d), v(v), rhs(rhs), area(area), p(p) {};

void Compartment::addChild(Compartment * child)
{
    branches.push_back(child);
};

void Compartment::addMechanismInstance(int mechType,
                                       double *  data, int dataSize,
                                       Datum  * pdata, int pdataSize)
{
    assert(mechanismsMap[mechType]!=-1);
    mechsTypes.push_back(mechType);
    if (dataSize>0)
        for (int i=0; i<dataSize; i++)
            this->data.push_back((floble_t) data[i]);
    if (pdataSize>0)
        for (int i=0; i<pdataSize; i++)
            this->pdata.push_back((offset_t) pdata[i]);
};

void Compartment::addVecPlay(double * t, double *y, PointProcInfo & ppi)
{
    assert(ppi.size>0);
    this->vecPlayInfo.push_back(ppi);
    for (int i=0; i<ppi.size; i++)
    {
        this->vecPlayTdata.push_back((floble_t) t[i]);
        this->vecPlayYdata.push_back((floble_t) y[i]);
    }
}

void Compartment::addNetCon(int preSynNrnThreadId, NetConX * nc, floble_t * weights)
{
    this->netconsPreSynIds.push_back(preSynNrnThreadId);
    this->netcons.push_back(NetConX(nc->mechType, nc->mechInstance, nc->delay,
                                   nc->weightIndex, nc->weightsCount, nc->active));
    for (int i=0; i<nc->weightsCount; i++)
        this->netconsWeights.push_back(weights[i]);

}
