#include "neurox/neurox.h"

#include <memory>
#include <utility>
#include <vector>

using namespace std;
using namespace neurox;
using namespace neurox::input;

class neurox::NetConX;

Compartment::~Compartment(){};

Compartment::Compartment(offset_t id, floble_t a, floble_t b, floble_t d,
                         floble_t v, floble_t rhs, floble_t area, offset_t p)
    : id(id), a(a), b(b), d(d), v(v), rhs(rhs), area(area), p(p){};

void Compartment::AddChild(Compartment *child) { branches.push_back(child); }

void Compartment::AddMechanismInstance(int mechType, int mechsInstance,
                                       double *data, int dataSize, Datum *pdata,
                                       int pdataSize) {
  assert(mechanisms_map[mechType] != -1);
  mechsTypes.push_back(mechType);
  mechsInstances.push_back(mechsInstance);
  if (dataSize > 0)
    for (int i = 0; i < dataSize; i++) this->data.push_back((floble_t)data[i]);
  if (pdataSize > 0)
    for (int i = 0; i < pdataSize; i++)
      this->pdata.push_back((offset_t)pdata[i]);
}

void Compartment::AddVecPlay(double *t, double *y, PointProcInfo &ppi) {
  assert(ppi.size > 0);
  this->vecPlayInfo.push_back(ppi);
  for (int i = 0; i < ppi.size; i++) {
    this->vecPlayTdata.push_back((floble_t)t[i]);
    this->vecPlayYdata.push_back((floble_t)y[i]);
  }
}

void Compartment::AddSerializedVdata(unsigned char *data, size_t size) {
  for (size_t i = 0; i < size; i++) this->vdata.push_back(data[i]);
}

void Compartment::AddNetCon(int preSynNrnThreadId, NetConX *nc,
                            floble_t *weights) {
  this->netconsPreSynIds.push_back(preSynNrnThreadId);
  this->netcons.push_back(NetConX(nc->mechType, nc->mechInstance, nc->delay,
                                  nc->weightIndex, nc->weightsCount,
                                  nc->active));
  for (int i = 0; i < nc->weightsCount; i++)
    this->netconsWeights.push_back(weights[i]);
}

void Compartment::ShrinkToFit() {
  branches.shrink_to_fit();
  mechsTypes.shrink_to_fit();
  mechsInstances.shrink_to_fit();
  data.shrink_to_fit();
  pdata.shrink_to_fit();
  vecPlayInfo.shrink_to_fit();
  vecPlayTdata.shrink_to_fit();
  vecPlayYdata.shrink_to_fit();
  vdata.shrink_to_fit();
  netcons.shrink_to_fit();
  netconsWeights.shrink_to_fit();
  netconsPreSynIds.shrink_to_fit();
}
