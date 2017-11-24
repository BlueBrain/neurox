#include "neurox/neurox.h"

#include <memory>
#include <utility>
#include <vector>

using namespace std;
using namespace neurox;
using namespace neurox::input;

class neurox::NetconX;

Compartment::~Compartment(){};

Compartment::Compartment(offset_t id, floble_t a, floble_t b, floble_t d,
                         floble_t v, floble_t rhs, floble_t area, offset_t p)
    : id_(id),
      a_(a),
      b_(b),
      d_(d),
      v_(v),
      rhs_(rhs),
      area_(area),
      p_(p),
      execution_time_(-1){};

void Compartment::AddChild(Compartment *child) { branches_.push_back(child); }

void Compartment::AddMechanismInstance(int mech_type, int mech_instances,
                                       double *data, int dataSize, Datum *pdata,
                                       int pdata_size) {
  assert(mechanisms_map_[mech_type] != -1);
  mechs_types_.push_back(mech_type);
  mechs_instances_.push_back(mech_instances);
  if (dataSize > 0)
    for (int i = 0; i < dataSize; i++) this->data.push_back((floble_t)data[i]);
  if (pdata_size > 0)
    for (int i = 0; i < pdata_size; i++)
      this->pdata.push_back((offset_t)pdata[i]);
}

void Compartment::AddVecplay(double *t, double *y, PointProcInfo &ppi) {
  assert(ppi.size > 0);
  this->vecplay_info_.push_back(ppi);
  for (int i = 0; i < ppi.size; i++) {
    this->vecplay_tdata_.push_back((floble_t)t[i]);
    this->vecplay_ydata_.push_back((floble_t)y[i]);
  }
}

void Compartment::AddSerializedVData(unsigned char *data, size_t size) {
  for (size_t i = 0; i < size; i++) this->vdata_.push_back(data[i]);
}

void Compartment::AddNetcon(int pre_syn_nrn_thread_id, NetconX *nc,
                            floble_t *weights) {
  this->netcons_pre_syn_ids_.push_back(pre_syn_nrn_thread_id);
  this->netcons_.push_back(NetconX(nc->mech_type_, nc->mech_instance_,
                                   nc->delay_, nc->weight_index_,
                                   nc->weights_count_, nc->active_));
  for (int i = 0; i < nc->weights_count_; i++)
    this->netcons_weights_.push_back(weights[i]);
}

void Compartment::ShrinkToFit() {
  branches_.shrink_to_fit();
  mechs_types_.shrink_to_fit();
  mechs_instances_.shrink_to_fit();
  data.shrink_to_fit();
  pdata.shrink_to_fit();
  vecplay_info_.shrink_to_fit();
  vecplay_tdata_.shrink_to_fit();
  vecplay_ydata_.shrink_to_fit();
  vdata_.shrink_to_fit();
  netcons_.shrink_to_fit();
  netcons_weights_.shrink_to_fit();
  netcons_pre_syn_ids_.shrink_to_fit();
}
