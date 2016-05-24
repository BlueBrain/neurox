#pragma once

#include "neurox/neurox.h"

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

    void addMechanism(int mechId, double * data, int dataSize, Datum * pdata, int pdataSize);

    vector<Compartment*> children;
    double a,b,d,v,rhs,area;
    vector<int> mechsIds;
    vector<double> data;
    vector<Datum> pdata;
  private:
};
