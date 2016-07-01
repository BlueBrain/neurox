#pragma once

#include "neurox/Neurox.h"

#include <memory>
#include <vector>

using namespace std;
using namespace Neurox;
using namespace Neurox::Input::Coreneuron;

Compartment::~Compartment()
{
    for (int c=0; c<branches.size(); c++)
        delete branches[c];
}

void Compartment::setValues(int id, double a, double b, double d, double v, double rhs, double area)
{
    this->id=id; this->a=a; this->b=b; this->d=d; this->v=v; this->rhs=rhs; this->area=area;
};

void Compartment::addChild(Compartment * child)
{
    branches.push_back(child);
};


void Compartment::addMechanism(int mechId, int instance, double * data, int dataSize, Datum * pdata, int pdataSize)
{
    mechsIds.push_back(mechId);
    mechsInstances.push_back(instance);
    this->data.insert (this->data.end() , &data[0] , &data[dataSize]);
    this->pdata.insert(this->pdata.end(), &pdata[0], &pdata[pdataSize]);
};
