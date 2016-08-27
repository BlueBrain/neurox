#pragma once

#include "neurox/neurox.h"
#include <memory>
#include <vector>

using namespace std;

namespace neurox
{
namespace Input
{
namespace Coreneuron
{
class Compartment
{
  public:
    Compartment()=delete;
    Compartment(int id, double a, double b, double d, double v, double rhs, double area, int p);
    ~Compartment();

    void addChild(Compartment* child);

    void addMechanismInstance(int mechId, double * data, int dataSize, Datum * pdata, int pdataSize);
    void addVecPlay(double * t, double *y, PointProcInfo & ppi);

    int id;
    vector<Compartment*> branches;
    double a,b,d,v,rhs,area;
    int p;
    vector<int> mechsTypes;
    vector<double> data;
    vector<Datum> pdata;

    //vecplay data
    vector<PointProcInfo> vecPlayInfo;
    vector<double> vecPlayTdata;
    vector<double> vecPlayYdata;

    //vdata
    vector<void*> vdata;
  private:
};

}; //Coreneuron
}; //Input
}; //NeuroX
