#pragma once

#include "neurox/Neurox.h"

#include "coreneuron/nrnconf.h"

#include <memory>
#include <vector>

using namespace std;

class Compartment
{
  public:
    Compartment()=delete;
    ~Compartment();

    void setSolverValues(double a, double b, double d, double v, double rhs, double area);

    void addChild(Compartment* child);

    void addMechanism(int mechId, int instance, double * data, int dataSize, Datum * pdata, int pdataSize);

    void addSynapse();

    vector<Compartment*> children;
    double a,b,d,v,rhs,area;
    vector<int> mechsIds;
    vector<int> mechsInstances;
    vector<double> data;
    vector<Datum> pdata;
  private:
};
