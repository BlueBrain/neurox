#pragma once

#include "neurox/Neurox.h"
#include <memory>
#include <vector>

using namespace std;

namespace Neurox
{

namespace Input
{

namespace Coreneuron
{

class Compartment
{
  public:
    Compartment()=delete;
    Compartment(int id, double a, double b, double d, double v, double rhs, double area);
    ~Compartment();

    void addChild(Compartment* child);

    void addMechanismInstance(int mechId, int instance, double * data, int dataSize, Datum * pdata, int pdataSize);

    int id;
    vector<Compartment*> branches;
    double a,b,d,v,rhs,area;
    vector<int> mechsIds;
    vector<int> mechsInstances;
    vector<double> data;
    vector<Datum> pdata;
    vector<NetConX> synapsesIn;
  private:
};

}; //Coreneuron
}; //Input
}; //Neurox
