/*
# =============================================================================
# Copyright (c) 2015 - 2021-2022 Blue Brain Project/EPFL
#
# See top-level LICENSE file for details.
# =============================================================================.
*/
#pragma once

#include "neurox/neurox.h"

#include <memory>
#include <vector>

using namespace std;

namespace neurox {
namespace input {
class Compartment {
 public:
  Compartment() = delete;
  Compartment(offset_t id, floble_t a, floble_t b, floble_t d, floble_t v,
              floble_t rhs, floble_t area, offset_t p);
  ~Compartment();

  void AddChild(Compartment* child);
  void AddMechanismInstance(int mech_type, int mech_instances, double* data,
                            int dataSize, Datum* pdata, int pdata_size);
  void AddVecplay(double* t, double* y, PointProcInfo& ppi);
  void AddSerializedVData(unsigned char* data, size_t size);
  void AddNetcon(int preSynNrnThreadId, NetconX* nc, floble_t* weights);
  void ShrinkToFit();

  offset_t id_;
  vector<Compartment*> branches_;
  floble_t a_, b_, d_, v_, rhs_, area_;
  offset_t p_;
  vector<int> mechs_types_;
  vector<int> mechs_instances_;
  vector<floble_t> data;
  vector<Datum> pdata;

  // vecplay data
  vector<PointProcInfo> vecplay_info_;
  vector<floble_t> vecplay_tdata_;
  vector<floble_t> vecplay_ydata_;

  // vdata (serialized)
  vector<unsigned char> vdata_serialized_;

  // netcons
  vector<NetconX> netcons_;
  vector<floble_t> netcons_weights_;
  vector<neuron_id_t> netcons_pre_syn_ids_;

 private:
};

};  // namespace input
};  // namespace neurox
