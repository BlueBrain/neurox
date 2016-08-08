#pragma once

#include "neurox/Neurox.h"

#include <memory>
#include <vector>

using namespace std;
using namespace NeuroX;
using namespace NeuroX::Input::Coreneuron;

Compartment::~Compartment(){};

Compartment::Compartment(int id, double a, double b, double d,
                            double v, double rhs, double area, int p):
    id(id), a(a), b(b), d(d), v(v), rhs(rhs), area(area), p(p) {};

void Compartment::addChild(Compartment * child)
{
    branches.push_back(child);
};


void Compartment::addMechanismInstance(int mechType, double * data, int dataSize, Datum * pdata, int pdataSize)
{
    mechsTypes.push_back(mechType);
    if (dataSize>0)
       this->data.insert (this->data.end() , data , data  + dataSize );
    if (pdataSize>0)
       this->pdata.insert (this->pdata.end() , pdata , pdata  + pdataSize );
};
