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
    Compartment(offset_t id, floble_t a, floble_t b, floble_t d,
                floble_t v, floble_t rhs, floble_t area, offset_t p);
    ~Compartment();

    void addChild(Compartment* child);

    void addMechanismInstance(int mechId,
                              double * data, int dataSize,
                              Datum * pdata, int pdataSize);
    void addVecPlay(double * t, double *y, PointProcInfo & ppi);

    offset_t id;
    vector<Compartment*> branches;
    floble_t a,b,d,v,rhs,area;
    offset_t p;
    vector<int> mechsTypes;
    vector<floble_t> data;
    vector<Datum> pdata;

    //vecplay data
    vector<PointProcInfo> vecPlayInfo;
    vector<floble_t> vecPlayTdata;
    vector<floble_t> vecPlayYdata;

    //vdata
    vector<void*> vdata;
  private:
};

}; //Coreneuron
}; //Input
}; //NeuroX
