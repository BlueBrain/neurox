#pragma once

#include "neurox/neurox.h"
#include <memory>
#include <vector>

using namespace std;

namespace neurox
{
namespace input
{
class Compartment
{
  public:
    Compartment()=delete;
    Compartment(offset_t id, floble_t a, floble_t b, floble_t d,
                floble_t v, floble_t rhs, floble_t area, offset_t p);
    ~Compartment();

    void AddChild(Compartment* child);
    void AddMechanismInstance(int mechId, int mechsInstance,
                              double * data, int dataSize,
                              Datum * pdata, int pdataSize);
    void AddVecPlay(double * t, double *y, PointProcInfo & ppi);
    void AddSerializedVdata(unsigned char * data, size_t size);
    void AddNetCon(int preSynNrnThreadId, NetConX * nc, floble_t * weights);
    void ShrinkToFit();

    offset_t id;
    vector<Compartment*> branches;
    floble_t a,b,d,v,rhs,area;
    offset_t p;
    vector<int> mechsTypes;
    vector<int> mechsInstances;
    vector<floble_t> data;
    vector<Datum> pdata;

    //vecplay data
    vector<PointProcInfo> vecPlayInfo;
    vector<floble_t> vecPlayTdata;
    vector<floble_t> vecPlayYdata;

    //vdata (serialized)
    vector<unsigned char> vdata;

    //netcons
    vector<NetConX > netcons ;
    vector<floble_t> netconsWeights;
    vector<neuron_id_t> netconsPreSynIds;
  private:
};

}; //Input
}; //NeuroX
