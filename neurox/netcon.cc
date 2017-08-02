#include "neurox/neurox.h"

#include <cstring>

using namespace neurox;

NetconX::NetconX() {}

NetconX::NetconX(int mech_type, offset_t mech_instance, floble_t delay,
                 offset_t weight_index, unsigned short weights_count,
                 bool active)
    : mech_type_(mech_type),
      weights_count_(weights_count),
      mech_instance_(mech_instance),
      delay_(delay),
      weight_index_(weight_index),
      active_(active) {}

NetconX::~NetconX() {}

void NetconX::Deliver(floble_t tt, Branch* branch) {
  if (this->active_)
    GetMechanismFromType(mech_type_)
        ->CallModFunction(branch, Mechanism::ModFunctions::kNetReceive, this,
                          tt);
}
