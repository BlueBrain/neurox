#pragma once

#include "neurox/Neurox.h"
#include <memory>
#include <vector>

using namespace std;

namespace NeuroX
{

namespace Input
{

namespace Coreneuron
{

class Compartment
{
  public:
    Compartment()=delete;
    Compartment(int id, double a, double b, double d, double v, double rhs, double area, int p=-999);
    ~Compartment();

    void addChild(Compartment* child);

    void addMechanismInstance(int mechId, double * data, int dataSize, Datum * pdata, int pdataSize, int pdataGap=0);

    int id;
    vector<Compartment*> branches;
    double a,b,d,v,rhs,area;
    int p;
    vector<int> mechsTypes;
    vector<double> data;
    vector<Datum> pdata;
    vector<NetConX> synapsesIn;
  private:
};

}; //Coreneuron
}; //Input
}; //NeuroX
