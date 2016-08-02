#pragma once

#include "neurox/Neurox.h"

#include <memory>
#include <vector>

using namespace std;
using namespace Neurox;
using namespace Neurox::Input::Coreneuron;

Compartment::~Compartment(){};

Compartment::Compartment(int id, double a, double b, double d,
                            double v, double rhs, double area):
    id(id), a(a), b(b), d(d), v(v), rhs(rhs), area(area) {};

void Compartment::addChild(Compartment * child)
{
    branches.push_back(child);
};


void Compartment::addMechanismInstance(int mechType, int instance, double * data, int dataSize, Datum * pdata, int pdataSize)
{
    mechsIds.push_back(mechType);
    mechsInstances.push_back(instance);

    if (data!=nullptr)
       this->data.insert (this->data.end() , data , data  + dataSize );

    if (pdata!=nullptr)
       this->pdata.insert(this->pdata.end(), pdata, pdata + pdataSize);
};
